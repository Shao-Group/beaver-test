#!/bin/bash
dir=`pwd`
STAR=$dir/../programs/STAR

# step 0: check all to-be-used tools/data
if [ "A" == "A" ];then
	echo "================================================================="
	echo "Check if to-be-used tools/datasets are properly linked..."
	echo "================================================================="
	if [ -e $STAR ];then
		echo -e "Tool STAR found successfully!"
	else
		echo -e "Tool STAR not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
		echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
    		exit 1
	fi
	echo -e "To-be-used tools/datasets found successfully!"
fi

# step 1: download reference
if [ "A" == "A" ];then
	cd $dir
	mkdir -p human
	cd human
	wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
	wget http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

	gzip -d Homo_sapiens.GRCh38.109.gtf.gz
	gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	gzip -d Homo_sapiens.GRCh38.cdna.all.fa.gz

	cd $dir
	mkdir -p mouse
	cd mouse
	wget http://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
	wget http://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
	wget http://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

	gzip -d Mus_musculus.GRCm39.109.gtf.gz
	gzip -d Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
	gzip -d Mus_musculus.GRCm39.cdna.all.fa.gz
fi

# step 2: indexing by STAR
if [ "A" == "A" ];then
	cd $dir/human
	mkdir -p human_index
	$STAR --runMode genomeGenerate \
     		--runThreadN 30 \
     		--genomeDir ./human_index \
     		--genomeFastaFiles $dir/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     		--sjdbGTFfile $dir/human/Homo_sapiens.GRCh38.109.gtf \
     		--sjdbOverhang 149

	cd $dir/mouse
	mkdir -p mouse_index
	$STAR --runMode genomeGenerate \
     		--runThreadN 30 \
     		--genomeDir ./mouse_index \
     		--genomeFastaFiles $dir/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa \
     		--sjdbGTFfile $dir/mouse/Mus_musculus.GRCm39.109.gtf \
     		--sjdbOverhang 149
fi

# indexing by Salmon, optional for additional experiments
if [ "A" == "B" ];then
	salmon=$dir/../programs/salmon
	cd $dir/human
	$salmon index -t Homo_sapiens.GRCh38.cdna.all.fa -i salmon.index -p 16
	cd $dir/mouse
	$salmon index -t Mus_musculus.GRCm39.cdna.all.fa -i salmon.index -p 16
fi
