
import sys
import numpy as np
from snputils import SNPObject
from snputils import VCFReader
from snputils import VCFWriter

def _chrom_to_int(chrom):
    chrom = chrom.replace('chr', '').replace('chrom', '')
    if chrom.isdigit():
        return int(chrom)
    chrom_map = {'X': 23, 'Y': 24, 'MT': 25}
    return chrom_map.get(chrom, 26)  # Assigns 26 to any other non-numeric chromosome


def convert_chromosome_X(snpobj: SNPObject) -> SNPObject:
    """ Transform chromosome X data.
    1. Convert all chromosome X identifiers to '23'.
    2. For male samples that only have haploid sequences on chromosome X, this function duplicates these sequences
    to simulate a diploid representation suitable for downstream analysis that expects diploid data.
    
    Args:
        snpobj (SNPObject): A SNPObject instance.
    
    Returns:
        SNPObject: A SNPObject with modified chromosome X data.
    """
    # Standardize chromosome labels for chromosome X to 23
    snpobj['variants_chrom'] = np.array([_chrom_to_int(chrom) for chrom in snpobj['variants_chrom']])
    
    maternal = snpobj.calldata_gt[:, :, 0]
    paternal = snpobj.calldata_gt[:, :, 1]

    # Chromosome X (23) variants
    chrom_X_idx = np.where(snpobj['variants_chrom'] == 23)[0]

    # Check samples where all values in the paternal strand are '-1' for chromosome X
    haploids = [col for col in range(paternal.shape[1]) if np.all(paternal[chrom_X_idx, col] == -1)]

    print(f"Samples identified for duplication in paternal data: {len(haploids)}")

    # Duplicate maternal data for haploid sequences in males
    for j in haploids:
        paternal[chrom_X_idx, j] = maternal[chrom_X_idx, j]

    # Update the paternal strand with the duplicated data
    snpobj.calldata_gt[:, :, 1] = paternal
    
    print(f"Processed {len(chrom_X_idx)} variants on chromosome X for {paternal.shape[1]} samples, duplicating haploid data in {len(haploids)} samples.")
    return snpobj


def main(input_vcf_path: str, output_vcf_path: str):
    """ Main function to process a VCF file, specifically targeting chromosome X data to adjust male haploid sequences.
    
    Args:
        input_vcf_path (str): Path to the input VCF file.
        output_vcf_path (str): Path to the output VCF file where the processed data will be saved.
    """
    print(f"Reading VCF data from {input_vcf_path}")
    reader = VCFReader(input_vcf_path)
    vcf = reader.read()

    print("Processing Chromosome X...")
    new_vcf = convert_chromosome_X(vcf)

    print(f"Writing processed VCF to {output_vcf_path}")
    vcf_writer = VCFWriter(new_vcf, output_vcf_path, phased=True)
    vcf_writer.write()
    print("Processing complete.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_chromosome_X.py <input_vcf_path> <output_vcf_path>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    main(input_vcf, output_vcf)
