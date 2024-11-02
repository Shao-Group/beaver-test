#!/bin/bash
dir=`pwd`
zumis=$dir/../programs/zUMIs/zUMIs.sh

# step 0: check all to-be-used tools/data
if [ "A" == "A" ];then
        echo "================================================================="
        echo "Check if to-be-used tools/datasets are properly linked..."
        echo "================================================================="
        if [ -e $zumis ];then
                echo -e "Tool zUMIs found successfully!"
        else
                echo -e "Tool zUMIs not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi
        echo -e "To-be-used tools/datasets found successfully!"
fi

# step 1: download
if [ "A" == "B" ];then
        cd $dir

	wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.R1.fastq.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.R2.fastq.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.I1.fastq.gz
	wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.Fibroblasts.GelCut.I2.fastq.gz
fi

# step 2: preprocess by tool zUMIs
if [ "A" == "A" ];then
        cd $dir
	sed -i "s! curDir! $dir! g" Mouse-Fibroblast.yaml
        $zumis -c -y Mouse-Fibroblast.yaml
fi

