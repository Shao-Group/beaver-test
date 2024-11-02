# Overview

This repository is dedicated to comparing the performance of the meta-assembler [**Aletsch**](https://github.com/Shao-Group/aletsch) against two leading meta-assemblers, [TransMeta](https://github.com/yutingsdu/TransMeta), [PsiCLASS](https://github.com/splicebox/PsiCLASS), as well as two meta-assembly pipelines based on single-sample assemblers, [StringTie2-Merge](https://ccb.jhu.edu/software/stringtie/index.shtml) and [Scallop2](https://github.com/Shao-Group/scallop2) combined with [TACO](https://tacorna.github.io). Here we provide instructions for downloading necessary tools, preparing datasets, executing the tools/pipelines, scoring Aletsch's output transcripts, and reproducing the results presented in the Aletsch paper.

# Step 1: Download and Link Tools

Our experiments involve the following tools:

| Tool                                                         | Version |                Description                |
| ------------------------------------------------------------ | :-----: | :---------------------------------------: |
| [Aletsch](https://github.com/Shao-Group/aletsch)             | v1.1.0  |              Meta Assembler               |
| [Transmeta](https://github.com/yutingsdu/TransMeta)          |  v.1.0  |              Meta Assembler               |
| [PsiCLASS](https://github.com/splicebox/PsiCLASS)            | v1.0.3  |              Meta Assembler               |
| [StringTie2](https://ccb.jhu.edu/software/stringtie/index.shtml) | v2.2.1  | Single-sample Assembler |
| [Scallop2](https://github.com/Shao-Group/scallop2)           | v1.1.2  |          Single-sample Assembler          |             |
| [STAR](https://github.com/alexdobin/STAR/tree/master)        | v2.7.11 |              RNA-seq Aligner              |
| [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#gffcompare_dl) | v0.11.2 |      Evaluate assembled transcripts       |
| [gtfcuff](https://github.com/Kingsford-Group/rnaseqtools)    |    -    |               RNA-seq tool                |

#### Step 1.1: Download Tools

* Access the homepages of the respective tools using the links provided above.

- Follow the download and compilation instructions on each tool's homepage.

#### Step 1.2: Link or Copy Executables

- For tools with available executable files, link or copy them to the `programs` directory. This includes `aletsch`, `scallop2`, `stringtie`, `STAR`, `gffcompare`, and `gtfcuff`.
- For tools without standalone executables (TransMeta and PsiCLASS), link the entire directory to `programs`.

Ensure the tools are accessible via the following paths:

```
your/path/to/programs/aletsch
your/path/to/programs/TransMeta/TransMeta
your/path/to/programs/psiclass/psiclass
your/path/to/programs/stringtie
your/path/to/programs/scallop2
your/path/to/programs/STAR
your/path/to/programs/gffcomapre
your/path/to/programs/gtfcuff
```

You may need to rename some of the executable files to match the paths listed above.

# Step 2: Download Datasets and Align

We evaluate the performance of the five methods using eight datasets, as outlined below. Each dataset is identified by its unique prefix (used in this repository) and accession ID for reference.

| Name in paper |      prefix in test(Ensembl)      |      Protocol       |                         Accession ID                         |
| :-----------: | :-------------------------------: | :-----------------: | :----------------------------------------------------------: |
|     SC-H1     |  **smartseq3_ensembl_human**  |      Smartseq3      | Random 100 cells from HEK293T of [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |            |
|     SC-M1     |    **smartseq3_ensembl_mouse**    |      Smartseq3      | All 369 cells from Mouse-Fibroblast of [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |

Use [STAR](https://github.com/alexdobin/STAR/tree/master) for read alignments for each sample/cell. For every dataset, compile a list of all BAM file paths as required by the different meta-assemblers. Example for Aletsch: `data/encode10_ensembl.star.list`. Your lists should follow this format: