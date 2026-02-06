"""
Design qPCR primers for genes from GTF annotations and genome FASTA.
"""

import argparse
import sys
from pathlib import Path

# Basics
import polars as pl

# Bioinformatics
import primer3
from Bio.Seq import Seq
from pyfaidx import Fasta

def import_gtf(fpath):
    # Import gene annotations
    gene_annot = pl.read_csv(fpath,
                             separator='\t', has_header=False)
    gene_annot = gene_annot.with_columns(pl.col('column_9').str.split('"').list.get(1).alias('column_10'))
    
    # Reindex start to zero
    gene_annot = gene_annot.with_columns((pl.col('column_4')-1))

    return gene_annot

def join_seq(df,
             chrom,
             fasta):
    ls = df.with_columns(
        pl.struct(['column_4', 'column_5']
                 ).map_elements(lambda x: fasta[chrom][x['column_4']:x['column_5']].seq
                               ).alias('sequence')
    )['sequence'].to_list()

    return ''.join(ls)

def pull_cds(seq_id,
             gene_annot,
             fasta,
             cds_only=False):
    '''
    seq_id : str
        Sequence ID of gene of interest.
    gene_annot : Polars DataFrame
        Dataframe of GTF file.
    fasta : pyfaidx object
        
    '''
    subdf = gene_annot.filter(pl.col('column_10') == seq_id)
    strand = list(subdf['column_7'])
    if not all([s==strand[0] for s in strand]):
        raise ValueError('Mixed strandedness for ID.')
    else:
        strand = strand[0]
        chrom = str(subdf['column_1'][0])

    if cds_only:
        full_mrna = subdf.filter(pl.col('column_3').is_in(["CDS"]))
    else:
        full_mrna = subdf.filter(pl.col('column_3').is_in(["5'-UTR","CDS","3'-UTR"]))
        
    full_seq = join_seq(full_mrna, chrom, fasta)

    if strand == '-':
        full_seq = str(Seq(full_seq).reverse_complement())
    else:
        pass

    return full_seq

def pick_primers(df,
                 pick,
                 allow_single=False,
                 name='pair'):

    # Set up coords to prevent overlap
    l_coords = []
    r_coords = []

    # Always take top primer pick
    pick_df = df[0]
    # Add coordinates of left primer to l_coords
    l_pos = pick_df['L-coords'][0][0]
    l_len = pick_df['L-coords'][0][1]
    l_coords.append(range(l_pos,l_pos+l_len))
    # Add coordinates of right primer to r_coords
    r_pos = pick_df['R-coords'][0][0]
    r_len = pick_df['R-coords'][0][1]
    r_coords.append(range(r_pos-r_len,r_pos+1))

    # Set a tracker for the number of primers
    tracker = 1
    for i in range(1,len(df)):
        # Check for overlapping intervals
        left = any([df[i]['L-coords'][0][0] in r for r in l_coords])
        right = any([df[i]['R-coords'][0][0] in r for r in r_coords])
        
        # Never allow both primers to be in the same interval
        if all([left,right]):
            pass
        # If one primer can be in the interval, then previous statement
        # has handled issues. Append.
        elif allow_single:
            pick_df = pl.concat([pick_df,df[i]])
            tracker += 1
        # If both primers must be unique, pass if either is in a previous
            # interval.
        elif any([left,right]):
            pass
        # If neither primer is in a previous interval, append.
        else:
            pick_df = pl.concat([pick_df, df[i]])
            tracker += 1

        # Check if the desired pair number has been hit.
        if tracker < pick:
            pass
        else:
            break

    # Check that the desired number of primers was found.
    if len(pick_df) == pick:
        pass
    else:
        print(f'Could not find {pick} primers in specification for {name}. Returning all {len(df)} primer pairs.')
        pick_df = df

    return pick_df
    

def design_primers(seq,
                   name='pair',
                   filt_gquad=True,
                   pick=2,
                   allow_single=False,):
    output = primer3.bindings.design_primers(seq_args={'SEQUENCE_TEMPLATE': (seq)},
                                global_args={
                                'PRIMER_OPT_SIZE': 20,
                                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                                'PRIMER_INTERNAL_MAX_SELF_END': 8,
                                'PRIMER_MIN_SIZE': 18,
                                'PRIMER_MAX_SIZE': 25,
                                'PRIMER_OPT_TM': 60.0,
                                'PRIMER_MIN_TM': 57.0,
                                'PRIMER_MAX_TM': 63.0,
                                'PRIMER_MIN_GC': 20.0,
                                'PRIMER_MAX_GC': 80.0,
                                'PRIMER_MAX_POLY_X': 100,
                                'PRIMER_INTERNAL_MAX_POLY_X': 100,
                                'PRIMER_SALT_MONOVALENT': 50.0,
                                'PRIMER_DNA_CONC': 50.0,
                                'PRIMER_MAX_NS_ACCEPTED': 0,
                                'PRIMER_MAX_SELF_ANY': 12,
                                'PRIMER_MAX_SELF_END': 8,
                                'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                                'PRIMER_PAIR_MAX_COMPL_END': 8,
                                'PRIMER_PRODUCT_SIZE_RANGE': [[100,220]]})

    left = output['PRIMER_LEFT']
    right = output['PRIMER_RIGHT']

    cols = ['PENALTY','SEQUENCE','COORDS','TM']

    df_dict = {'L-sequence':[],
               'L-penalty':[],
               'L-tm':[],
               'L-coords':[],
               'R-sequence':[],
               'R-penalty':[],
               'R-tm':[],
               'R-coords':[]}
    for i,v in enumerate(left):
        r_dict = right[i]
        for c in cols:
            df_dict[f'L-{c.lower()}'].append(v[c])
            df_dict[f'R-{c.lower()}'].append(r_dict[c])

    df = pl.DataFrame(df_dict)
    if filt_gquad:
        df = df.filter(~pl.col('L-sequence').str.contains('GGGG'))
        df = df.filter(~pl.col('R-sequence').str.contains('GGGG'))

    if not pick:
        pass
    elif len(df) < pick:
        print(f'Fewer than {pick} total primers for {name}. Returning all primers.')
    else:
        df = pick_primers(df,pick,allow_single,name)
        
    # Name output and get amplicon size
    df = df.with_columns(
        pl.format(f'{name}_{{}}', pl.int_range(pl.len())).alias('name'))
    df = df.with_columns((pl.col('R-coords').list.get(0) -
                          pl.col('L-coords').list.get(0)
                         ).alias('amplicon-bp'))
    df = df.select(['name','amplicon-bp']+list(df_dict.keys()))
    
    return df

def full_qpcr(gene_id,
              gene_annot,
              fasta,
              cds_only=False,
              filt_gquad=True,
              pick=2,
              allow_single=False,):
    '''
    Pull CDS for gene and design primers

    Parameters:
    -----------
    gene_id : str
        Name of transcript, as named in GTF file.
    gene_annot : Polars DataFrame
        DataFrame of GTF file. Required columns:
        - "column_1"  : chromosome, name matching FASTA
        - "column_3"  : gene feature
        - "column_4"  : start coordinate
        - "column_5"  : end coordinate
        - "column_7"  : strand
        - "column_10" : transcript ID
    fasta : pyfaidx object
        Full genome, imported using pyfaidx (required for fast indexing)
    cds_only : Bool
        If True, only evaluates CDS. If False, also includes UTRs. Default: False
    filt_gquad : Bool
        If True, removes any primers with G-quadruplex (GGGG). Default: True
    pick : int
        Number of primers to return after filtering. Default: 2
    allow_single : Bool
        If True, returned primer pairs may share either the forward or reverse primer
        coordinates (but not both). If False, no primer pairs may share forward or 
        reverse coordinates.

    Returns:
    --------
    A dataframe containing qPCR primers with amplicons 100-220 bp and Tm ~60°C.
    '''
    seq = pull_cds(gene_id,
                   gene_annot,
                   fasta,
                   cds_only)
    
    primers = design_primers(seq,
                             name=gene_id,
                             filt_gquad=filt_gquad,
                             pick=pick,
                             allow_single=allow_single,)

    return primers

def main():
    parser = argparse.ArgumentParser(
        description='Design qPCR primers for specified genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -g annotations.gtf -f genome.fa -i gene1 gene2 gene3
  %(prog)s -g annotations.gtf -f genome.fa -i gene1 --cds_only --pick 5
  %(prog)s -g annotations.gtf -f genome.fa -i gene1 gene2 --output primers.csv
        """
    )
    
    # Required arguments
    parser.add_argument('-g', '--gtf', 
                        required=True,
                        type=str,
                        help='Path to GTF annotation file')
    
    parser.add_argument('-f', '--fasta',
                        required=True, 
                        type=str,
                        help='Path to genome FASTA file')
    
    parser.add_argument('-i', '--genes',
                        required=True,
                        nargs='+',
                        type=str,
                        help='Gene/transcript IDs to design primers for (space-separated)')
    
    # Optional arguments
    parser.add_argument('-o', '--output',
                        type=str,
                        default=None,
                        help='Output CSV file (default: print to stdout)')
    
    parser.add_argument('--cds_only',
                        action='store_true',
                        default=False,
                        help='Use only CDS regions (default: include UTRs)')
    
    parser.add_argument('--filt_gquad',
                        action='store_true',
                        default=True,
                        help='Filter out primers with G-quadruplex motifs (default: True)')
    
    parser.add_argument('--no_filt_gquad',
                        dest='filt_gquad',
                        action='store_false',
                        help='Do not filter G-quadruplex motifs')
    
    parser.add_argument('--pick',
                        type=int,
                        default=2,
                        help='Number of primer pairs to return (default: 2)')
    
    parser.add_argument('--allow_single',
                        action='store_true',
                        default=False,
                        help='Allow primer pairs to share forward or reverse coordinates')
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.gtf).exists():
        sys.stderr.write(f"Error: GTF file not found: {args.gtf}\n")
        sys.exit(1)
    
    if not Path(args.fasta).exists():
        sys.stderr.write(f"Error: FASTA file not found: {args.fasta}\n")
        sys.exit(1)
    
    # Load annotations and genome
    print(f"Loading GTF annotations from {args.gtf}...", file=sys.stderr)
    gene_annot = import_gtf(args.gtf)
    
    print(f"Loading genome FASTA from {args.fasta}...", file=sys.stderr)
    fasta = Fasta(args.fasta)
    
    # Design primers for each gene
    all_primers = []
    for gene_id in args.genes:
        print(f"Designing primers for {gene_id}...", file=sys.stderr)
        try:
            primers = full_qpcr(
                gene_id,
                gene_annot,
                fasta,
                cds_only=args.cds_only,
                filt_gquad=args.filt_gquad,
                pick=args.pick,
                allow_single=args.allow_single
            )
            all_primers.append(primers)
        except Exception as e:
            print(f"Warning: Failed to design primers for {gene_id}: {e}", 
                  file=sys.stderr)
            continue
    
    # Combine results
    if not all_primers:
        sys.stderr.write("Error: No primers were successfully designed\n")
        sys.exit(1)
    
    result = pl.concat(all_primers)
    for c in result.columns:
        if 'coords' in c:
            result = result.with_columns(pl.col(c).list.get(0))
        elif 'tm' in c:
            result = result.with_columns(pl.col(c).round(1))
        elif 'penalty' in c:
            result = result.with_columns(pl.col(c).round(4))
        else:
            pass
    
    # Output results
    if args.output:
        result.write_csv(args.output)
        print(f"Primers written to {args.output}", file=sys.stderr)
    else:
        print(result.write_csv())
    
    print(f"Successfully designed primers for {len(all_primers)} gene(s)", 
          file=sys.stderr)


if __name__ == '__main__':
    main()