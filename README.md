# Overview

This repository is dedicated to comparing the performance of our new assembler [**Beaver**](https://github.com/Shao-Group/beaver) against three leading meta-assemblers, [Aletsch](https://github.com/Shao-Group/aletsch), [TransMeta](https://github.com/yutingsdu/TransMeta), [PsiCLASS](https://github.com/splicebox/PsiCLASS), as well as two popular single-sample assemblers, [StringTie2](https://ccb.jhu.edu/software/stringtie/index.shtml) and [Scallop2](https://github.com/Shao-Group/scallop2). Here we provide instructions for downloading necessary tools, preparing datasets, executing the tools/pipelines, scoring Beaver's output transcripts, and reproducing the results presented in the Beaver paper.

# Step 1: Download and Link Tools

Our experiments involve the following tools:

| Tool                                                         | Version |                Description                |
| ------------------------------------------------------------ | :-----: | :---------------------------------------: |
| [Beaver](https://github.com/Shao-Group/beaver)             | v1.0.0  |              Cell-specific Assembler at single-cell resolution              |
| [Aletsch](https://github.com/Shao-Group/aletsch)             | v1.1.0  |              Meta Assembler               |
| [Transmeta](https://github.com/yutingsdu/TransMeta)          |  v.1.0  |              Meta Assembler               |
| [PsiCLASS](https://github.com/splicebox/PsiCLASS)            | v1.0.3  |              Meta Assembler               |
| [StringTie2](https://ccb.jhu.edu/software/stringtie/index.shtml) | v2.2.1  | Single-sample Assembler |
| [Scallop2](https://github.com/Shao-Group/scallop2)           | v1.1.2  |          Single-sample Assembler          |             |
| [STAR](https://github.com/alexdobin/STAR/tree/master)        | v2.7.11 |              RNA-seq Aligner              |
| [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#gffcompare_dl) | v0.11.2 |      Evaluate assembled transcripts       |

#### Step 1.1: Download Tools

* Access the homepages of the respective tools using the links provided above.

- Follow the download and compilation instructions on each tool's homepage.

#### Step 1.2: Link or Copy Executables

- For tools with available executable files, link or copy them to the `programs` directory. This includes `beaver`, `aletsch`, `scallop2`, `stringtie`, `STAR`, `gffcompare`, and `gtfcuff`.
- For tools without standalone executables (TransMeta and PsiCLASS), link the entire directory to `programs`.

Ensure the tools are accessible via the following paths:

```
your/path/to/programs/beaver
your/path/to/programs/aletsch
your/path/to/programs/TransMeta/TransMeta
your/path/to/programs/psiclass/psiclass
your/path/to/programs/stringtie
your/path/to/programs/scallop2
your/path/to/programs/STAR
your/path/to/programs/gffcomapre
```

You may need to rename some of the executable files to match the paths listed above.

# Step 2: Download Datasets and Align

We evaluate the performance of the six methods using four datasets, as outlined below. Each dataset is identified by its unique prefix (used in this repository) and accession ID for reference.

| Dataset |      # Cells     |      Protocol       |                         Accession ID                         |
| :-----------: | :-------------------------------: | :-----------------: | :----------------------------------------------------------: |
|     HEK293T     |  192  |      Smart-seq3      | [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |            |
|     Mouse-Fibroblast     |    369    |      Smart-seq3      | [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |

Use [STAR](https://github.com/alexdobin/STAR/tree/master) for read alignments for each sample/cell. For every dataset, compile a list of all BAM file paths as required by the different meta-assemblers. 