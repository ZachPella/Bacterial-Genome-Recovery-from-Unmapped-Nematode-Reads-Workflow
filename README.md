# Metagenome Assembly and Binning Pipeline

A comprehensive metagenomics pipeline for bacterial genome recovery from unmapped sequencing reads. This workflow performs *de novo* assembly, genome binning, quality assessment, taxonomic classification, and functional annotation of microbial genomes.

## Overview

This pipeline takes unmapped reads from host-associated sequencing data and reconstructs bacterial genomes through metagenome assembly and binning. The workflow follows best practices for metagenome-assembled genome (MAG) recovery and includes comprehensive quality control and annotation steps.

## Pipeline Steps

### 1. Extract Unmapped Reads (`1_get_unmapped_reads`)
- Extracts paired-end reads where both mates are unmapped from BAM files
- Converts to FASTQ format for metagenome assembly
- **Input**: `*.sorted.refrename.bam`
- **Output**: Unmapped FASTQ files

### 2. Taxonomic Classification (`2_kraken2`)
- Initial taxonomic profiling using Kraken2
- Identifies microbial diversity before assembly
- **Output**: Taxonomic composition reports

### 3. Extract Bacterial Reads (`3_extract_bacteria_only`)
- Filters reads classified as bacterial origin
- Reduces assembly complexity by removing host/viral/archaeal contamination
- **Output**: Bacterial-specific FASTQ files

### 4. Metagenome Assembly (`4_megahit_contigs`)
- *De novo* assembly using MEGAHIT
- Optimized for metagenome data with varying coverage
- **Output**: Assembled contigs in FASTA format

### 5. Binning Preparation (`5_binning_prep`)
- Maps reads back to assembled contigs
- Calculates coverage depth across contigs
- Prepares input files for genome binning
- **Output**: Sorted BAM files and coverage tables

### 6. Genome Binning (`6_metabat`)
- Bins contigs into putative genomes using MetaBAT2
- Uses tetranucleotide frequency and coverage patterns
- **Output**: Individual genome bins (FASTA files)

### 7. Quality Assessment (`7_checkm`)
- Evaluates genome completeness and contamination using CheckM
- Assesses bin quality based on single-copy marker genes
- **Output**: Quality metrics for each bin

### 8. Taxonomic Classification (`8_gtdbtk`)
- Assigns taxonomy to genome bins using GTDB-Tk
- Places genomes in standardized bacterial taxonomy
- **Output**: Taxonomic assignments and phylogenetic placement

### 9. Functional Annotation (`9_annotation`)
- Annotates protein-coding genes in high-quality bins
- Predicts gene function and metabolic pathways
- **Output**: Annotated genomes with functional predictions

## Requirements

### Software Dependencies
- **SLURM** (job scheduler)
- **samtools** (≥1.9)
- **Kraken2** with bacterial database
- **MEGAHIT** (≥1.2.9)
- **BWA-MEM** or **Bowtie2**
- **MetaBAT2** (≥2.15)
- **CheckM** (≥1.1.3)
- **GTDB-Tk** (≥2.0) with GTDB database
- **Prokka** or **Prodigal** + annotation tools

### Database Requirements
- **Kraken2 database** (standard or custom bacterial)
- **CheckM database** (marker gene sets)
- **GTDB-Tk database** (≥200GB, latest release)
- **Annotation databases** (COG, KEGG, Pfam, etc.)

### System Requirements
- HPC cluster with SLURM
- Memory: 16-200GB depending on step
- CPUs: 4-32 cores per job
- Storage: ~100GB per sample for intermediate files
- High-memory nodes recommended for GTDB-Tk

## Usage

### Sequential Workflow
```bash
# 1. Extract unmapped reads
sbatch 1_get_unmapped_reads

# 2. Initial taxonomic profiling
sbatch 2_kraken2

# 3. Extract bacterial reads only
sbatch 3_extract_bacteria_only

# 4. Metagenome assembly
sbatch 4_megahit_contigs

# 5. Prepare for binning
sbatch 5_binning_prep

# 6. Genome binning
sbatch 6_metabat

# 7. Quality assessment
sbatch 7_checkm

# 8. Taxonomic classification
sbatch 8_gtdbtk

# 9. Functional annotation
sbatch 9_annotation
```

### Input Data Structure
```
project_directory/
├── *.sorted.refrename.bam    # Input BAM files
├── 1_get_unmapped_reads      # Pipeline scripts
├── 2_kraken2
├── ...
└── 9_annotation
```

### Output Structure
```
project_directory/
├── unmapped_reads/          # Extracted unmapped reads
├── kraken2_results/         # Initial taxonomic profiles
├── bacterial_reads/         # Filtered bacterial reads
├── assembly/
│   ├── contigs.fa          # Assembled contigs
│   └── assembly_stats.txt  # Assembly metrics
├── mapping/                 # Read mapping results
├── bins/                    # Individual genome bins
├── checkm_results/          # Quality assessment
├── gtdbtk_results/         # Taxonomic assignments
└── annotation/             # Functional annotations
```

## Key Features

- **Comprehensive**: Complete workflow from reads to annotated genomes
- **Quality-focused**: Multiple QC steps ensure high-quality results
- **Scalable**: Optimized for HPC environments
- **Standardized**: Uses established tools and databases
- **Reproducible**: Clear documentation and version control

## Expected Results

### Assembly Metrics
- **N50**: Typically 1-50kb for metagenomes
- **Total length**: Varies by microbial diversity
- **Contig count**: 1,000-100,000+ depending on complexity

### Genome Bins
- **High-quality bins**: >90% complete, <5% contamination
- **Medium-quality bins**: >50% complete, <10% contamination
- **Typical yield**: 5-50 bins per sample depending on diversity

### Taxonomic Coverage
- **Novel taxa**: May recover previously uncultured bacteria
- **Known groups**: Confirm presence of expected bacterial families
- **Phylogenetic placement**: High-resolution taxonomic assignments

## Quality Thresholds

### CheckM Quality Standards
- **High-quality**: ≥90% complete, ≤5% contamination
- **Medium-quality**: ≥50% complete, ≤10% contamination
- **Low-quality**: <50% complete or >10% contamination

### Assembly Quality
- **Minimum contig length**: 1,000 bp
- **N50 target**: >10 kb for good assemblies
- **Coverage depth**: >10× for reliable binning

## Customization

### Adjusting Parameters
- **Assembly**: Modify k-mer sizes in MEGAHIT
- **Binning**: Adjust minimum contig length for MetaBAT2
- **Quality**: Change completeness/contamination thresholds
- **Annotation**: Select specific databases for functional annotation

### Database Updates
- **GTDB-Tk**: Update to latest release for current taxonomy
- **CheckM**: Ensure compatible marker gene sets
- **Functional**: Update COG, KEGG, Pfam databases

## Performance Optimization

### Memory Management
- **Assembly**: 32-64GB for complex metagenomes
- **GTDB-Tk**: 200GB+ RAM recommended
- **Binning**: 16-32GB typically sufficient

### Parallel Processing
- **Assembly**: Use all available cores
- **Mapping**: Parallelize across samples
- **Annotation**: Process bins in parallel

## Troubleshooting

### Common Issues
1. **Assembly failures**: Reduce k-mer size or increase memory
2. **Low bin quality**: Check coverage depth and assembly quality
3. **GTDB-Tk errors**: Verify database installation and paths
4. **Annotation timeouts**: Increase time limits for large genomes

### Quality Control Checkpoints
- **Assembly stats**: Check N50 and total length
- **Mapping rates**: Ensure >80% reads map back
- **Bin statistics**: Review completeness distributions
- **Taxonomy**: Verify reasonable phylogenetic placements

## Citation

If you use this pipeline, please cite the relevant tools:
- **MEGAHIT**: Li et al. (2015) Bioinformatics
- **MetaBAT2**: Kang et al. (2019) PeerJ
- **CheckM**: Parks et al. (2015) Genome Research
- **GTDB-Tk**: Chaumeil et al. (2020) Bioinformatics

## Best Practices

### Sample Preparation
- **Read depth**: >10M paired reads recommended
- **Quality**: Remove adapters and low-quality reads
- **Contamination**: Screen for host DNA removal

### Analysis Strategy
- **Co-assembly**: Consider pooling related samples
- **Validation**: Cross-reference results across tools
- **Documentation**: Track versions and parameters used

## License

MIT License - free to use and modify for research purposes.

---

**Pipeline for metagenome-assembled genome (MAG) recovery**  
*From unmapped reads → High-quality bacterial genomes*
