#!/bin/bash
#rsem=/datadisk1/qqz5133/tools/RSEM/software/bin
#PG=/datadisk1/qqz5133/nf-merge/link-test/programs
dir=`pwd`
rawdata=$dir/Mouse-Fibroblast_index
refgtf=/data/shared/data/mouse/Mus_musculus.GRCm39.110.gtf
#reffa=$dir/../data/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa
#simfa=$dir/Mouse-Fibroblast_simulation/fastq
#simml=$dir/Mouse-Fibroblast_simulation/model
#simr=$dir/Mouse-Fibroblast_simulation/simulated_fastq
simbam=$dir/Mouse-Fibroblast_simulation/bam
simgt=$dir/Mouse-Fibroblast_simulation/ground_truth
num=369
ref=$dir/mouse_ref/mouse_rsem
starref=$dir/mouse_ref

testlist='2 5 10 30 50 80 100 150 200 300 369'
#testlist='200 300 369'
result=$dir/Mouse-Fibroblast_simulation/results
transmeta=$dir/../programs/TransMeta/TransMeta
psiclass=$dir/../programs/psiclass/psiclass
gffcompare=$dir/../programs/gffcompare
gtfcuff=$dir/../programs/gtfcuff
lmversion=8
linkmerge=$dir/../programs/scAletsch-${lmversion}
scallop2=$dir/../programs/scallop2
stringtie2=$dir/../programs/stringtie2

# step 1: prepare reference
if [ "A" == "B" ];then
	mkdir -p $dir/mouse_ref
	$rsem/rsem-prepare-reference --gtf $refgtf $reffa $ref
fi

# step 1.1: bam to fastq
if [ "A" == "B" ];then
	mkdir -p $simfa

	joblist=$simfa/job.list
	rm -rf $joblist

	for ((i=1;i<=$num; i++))
        do
		echo "samtools fastq -1 $simfa/${i}_1.fastq -2 $simfa/${i}_2.fastq $rawdata/$i.bam" >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 16 bash -c CMD 1> /dev/null 2> /dev/null

fi

# step 2: learn real data
if [ "A" == "B" ];then
	mkdir -p $simml
	mkdir -p $simr

	joblist=$dir/job.list
	rm -rf $joblist

	for ((i=1;i<=$num; i++))
        do
		echo "$rsem/rsem-calculate-expression --star --star-path $PG --paired-end --estimate-rspd --append-names --output-genome-bam -p 2 --single-cell-prior --calc-pme $simfa/${i}_1.fastq $simfa/${i}_2.fastq $ref $simml/$i" >> $joblist
	done 
	
	cat $joblist | xargs -L 1 -I CMD -P 20 bash -c CMD 1> /dev/null 2> /dev/null
fi

# step 3: simulate 
if [ "A" == "B" ];then
	mkdir -p $simr
	joblist=$simr/job.list
	rm -rf $joblist

	for ((i=1;i<=$num; i++))
        do
		lines="$(wc -l $simfa/${i}_1.fastq | awk '{print $1}')"
		reads="$(echo $((lines / 4)))"

		#theta=$(sed '3q;d' $simml/$i.stat/$i.theta | awk '{print $1}')
		theta=0.1

		echo "$rsem/rsem-simulate-reads $ref $simml/$i.stat/$i.model $simml/$i.isoforms.results $theta $reads $simr/$i --seed 0" >> $joblist
	done

	cat $joblist | xargs -L 1 -I CMD -P 20 bash -c CMD 1> /dev/null 2> /dev/null
fi

# step 4: align to generate bam
if [ "A" == "B" ];then
        mkdir -p $simbam

	tnum=4

        joblist=$simbam/job.list
        rm -rf $joblist

        for ((i=1;i<=$num; i++))
        do
                script=$simbam/$i.sh
                echo "#!/bin/bash" > $script
                echo "mkdir -p $simbam/$i" >> $script
                echo "cd $simbam/$i" >> $script
                echo "${PG}/STAR --runThreadN $tnum --outSAMstrandField intronMotif --chimSegmentMin 20 --genomeDir $starref --twopassMode Basic --readFilesIn $simr/${i}_1.fq $simr/${i}_2.fq" >> $script
                echo "samtools view --threads $tnum -b Aligned.out.sam > Aligned.out.bam" >> $script
                echo "samtools sort --threads $tnum Aligned.out.bam > Aligned.out.sort.bam" >> $script
                echo "mv Aligned.out.sort.bam $simbam/$i.bam" >> $script

                chmod +x $script
                echo $script >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 5 bash -c CMD 1> /dev/null 2> /dev/null

fi

# step 5: get isoform ground truth
if [ "A" == "B" ];then
        mkdir -p $simgt

        joblist=$simr/gt.job.list
        rm -rf $joblist

        for ((i=1;i<=$num; i++))
        do
                ifile=$simr/$i.sim.isoforms.results
                ofile=$simgt/$i.list
                echo "python $dir/get_ground_truth.py $i $ifile $ofile" >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 20 bash -c CMD 1> /dev/null 2> /dev/null
fi

ref=$refgtf
index=$simbam
# step 6: run transmeta and psiclass
if [ "A" == "B" ];then
        echo "Running Transmeta..."

        cur=$result/transmeta
        mkdir -p $cur
        cd $cur

        # generate job list
        joblist=$cur/jobs.list
        rm -rf $joblist

        for num in $testlist;
        do
                list=$cur/$num.list
                rm -rf $list

                for ((i=1;i<=$num; i++))
                do
                        echo "$index/$i.bam" >> $list
                done

                mkdir -p $cur/merge_$num

                script=$cur/$num.sh
                rm -rf $script

                echo "$transmeta -B $list -s unstranded -o $cur/merge_$num -p 2 > $cur/merge_$num/merge_$num.log" > $script
                # gffcompare individual-level
                for ((j=1;j<=$num; j++))
                do
                        echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/individual.${j}.stats $cur/merge_$num/TransMeta.bam$j.gtf" >> $script
                        echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $cur/merge_$num/TransMeta.bam$j.gtf > $cur/merge_$num/TransMeta.bam$j.otherchrm.gtf" >> $script
                        echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/individual.${j}.other.stats $cur/merge_$num/TransMeta.bam$j.otherchrm.gtf" >> $script
                done

                chmod +x $script
                echo $script >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 10 bash -c CMD 1> /dev/null 2> /dev/null &

        echo "Finish Transmeta."
fi

# PsiCLASS
if [ "A" == "B" ];then
        echo "Running PsiCLASS..."

        cur=$result/psiclass
        mkdir -p $cur
        cd $cur

        # generate job list
        joblist=$cur/jobs.list
        rm -rf $joblist

        for num in $testlist;
        do
                list=$cur/$num.list
                rm -rf $list

                for ((i=1;i<=$num; i++))
                do
                        echo "$index/$i.bam" >> $list
                done

                mkdir -p $cur/merge_$num

                script=$cur/$num.sh
                rm -rf $script

                for cov in '1';
                do
                        mkdir -p $cur/merge_${num}/cov_$cov
                        echo "$psiclass --lb $list -o $cur/merge_$num/cov_$cov/ --vd $cov -p 2 > $cur/merge_$num/cov_$cov/merge_${num}_$cov.log" > $script
                        # gffcompare individual-level
                        for ((j=0;j<$num; j++))
                        do
                                echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/cov_1/individual.${j}.stats $cur/merge_$num/cov_1/psiclass_sample_$j.gtf" >> $script
                                echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $cur/merge_$num/cov_1/psiclass_sample_$j.gtf > $cur/merge_$num/cov_1/psiclass_sample_$j.otherchrm.gtf" >> $script
                                echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/individual.${j}.other.stats $cur/merge_$num/cov_1/psiclass_sample_$j.otherchrm.gtf" >> $script
                        done
                done

                chmod +x $script
                echo $script >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 10 bash -c CMD 1> /dev/null 2> /dev/null &

        echo "Finish PsiCLASS."
fi

# Aletsch
aletschExe=aletsch-585c8a 
aletsch=$dir/../programs/$aletschExe
#model=$dir/../models/aletsch-39.ensembl.chr1.joblib
model=$dir/../models/aletsch-rsemSim_mouse.joblib

if [ "A" == "B" ];then
        echo "Running Aletsch..."

        cur=$result/aletsch
        rm -rf $cur
        mkdir -p $cur
        cd $cur

        # generate job list
        joblist=$cur/jobs.list
        rm -rf $joblist

        # generate bai
        #cd $index
        #for ((i=1;i<=369; i++))
        #do
                #samtools index $index/$i.bam
                #mv $index/$i.bam.bai $index/$i.bai
        #done

        cd $cur
        for num in $testlist;
	do
                list=$cur/$num.list
                rm -rf $list

                for ((i=1;i<=$num; i++))
                do
                    echo "$index/$i.bam $index/$i.bai paired_end" >> $list
                done

                test=$cur/merge_$num
                rm -rf $test
                mkdir -p $test

                script=$cur/$num.sh
                rm -rf $script

                profile=$test/profile
                mkdir -p $profile

                sgtf=$test/gtf
                mkdir -p $sgtf
                rm -rf $sgtf/gff-scripts

                #echo "source ~/.profile" > $script
                echo "cd $test" > $script

                echo "$aletsch --profile -i $list -p profile > profile.preview" >> $script
                echo "{ /usr/bin/time -v $aletsch -i $list -o meta.gtf -d gtf -p $profile -c $((2 * $num)) > aletsch.log ; } 2> merge_$num.time" >> $script
                echo "$gffcompare -M -N -r $ref -o gffmul.stats meta.gtf" >> $script
                #echo "refnum=\`cat gffmul.stats | grep Reference | grep mRNA | head -n 1 | awk '{print \$5}'\`" >> $script
		#echo "$gtfcuff roc gffmul.meta.gtf.tmap \$refnum cov > roc.mul" >> $script
                #echo "/datadisk1/qzs23/tools/conda/env/qzs23/bin/python $dir/score_individual.py -i $sgtf -m $model -c $num -o $sgtf -p 0 > score.log" >> $script

                maxid=`cat $list | wc -l`
                ln -sf $ref $sgtf
                ln -sf $gffcompare $sgtf
                ln -sf $gtfcuff $sgtf
                for k in `seq 0 $maxid`
                do
                        kscript=$sgtf/$k.sh
                        echo "./gffcompare -M -N -r `basename $ref` -o $k.mul.stats $k.gtf" > $kscript
                        #echo "refnum=\`cat $k.mul.stats | grep Reference | grep mRNA | head -n 1 | awk '{print \$5}'\`" >> $kscript
                        #echo "./gtfcuff roc $k.mul.$k.gtf.tmap \$refnum cov > $k.roc.mul" >> $kscript
                        echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $k.gtf > $k.otherchrm.gtf" >> $kscript
                        echo "$gffcompare -r $ref -o $k.other.stats $k.otherchrm.gtf" >> $kscript
                        echo "$kscript" >> $sgtf/gff-scripts
                        chmod u+x $kscript
                done

                echo "cd $sgtf" >> $script
                echo "cat gff-scripts | xargs -L 1 -I CMD -P 10 bash -c CMD 1> /dev/null 2> /dev/null &" >> $script
                chmod +x $script
                echo $script >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 10 bash -c CMD 1> /dev/null 2> /dev/null &
	echo "Finish Aletsch."
fi

# Aletsch Scoring
if [ "A" == "B" ];then
        echo "Running Aletsch scoring..."

        cur=$result/aletsch
        cd $cur

        for num in $testlist;
        do
            test=$cur/merge_$num
            sgtf=$test/gtf
            #echo "python3 $dir/score_individual.py -i $sgtf -m $model -c $num -o $sgtf -p 0 > $test/score.log"
            nohup /home/grads/qzs23/.conda/envs/qzs23-env/bin/python $dir/score_individual.py -i $sgtf -m $model -c $num -o $sgtf -p 0 > $test/score.log 2>&1 &
        done
        
        echo "Finish Aletsch scoring."
fi

# Aletsch-link
if [ "A" == "A" ];then
        echo "Running altesch - Link..."

        sgtf=$result/aletsch-link-${lmversion}/raw-assemblies
        cur=$result/aletsch-link-${lmversion}
        mkdir -p $cur
        cd $cur

        #rm -rf $sgtf
        mkdir -p $sgtf

        # generate job list
        joblist=$cur/jobs.list
        rm -rf $joblist

        for num in $testlist;
        do
                list=$cur/$num.list
                rm -rf $list

                mkdir -p $cur/merge_$num
                mkdir -p $sgtf/merge_$num
                for ((i=1;i<=$num; i++))
                do
                        cp $result/aletsch/merge_$num/gtf/$(($i - 1)).prob.gtf $sgtf/merge_$num/$i.gtf
                        echo "$sgtf/merge_$num/$i.gtf" >> $list
                done

                script=$cur/$num.sh
                rm -rf $script
                echo "$linkmerge $list $cur/merge_$num/linkmerge_$num > $cur/merge_$num/merge_$num.log" > $script
                echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/merge_${num}.stats $cur/merge_$num/linkmerge_$num.gtf" >> $script
                #echo "$gtfcuff roc $cur/merge_$num/merge_${num}.linkmerge_$num.gtf.tmap 225036 cov > $cur/merge_$num/linkmerge_$num.roc.txt" >> $script
		
		# gffcompare individual-level
                lmsgtf=$cur/merge_$num/linkmerge_${num}_sgtf
                mkdir -p $lmsgtf
                cd $lmsgtf
                rm -rf gff-scripts
                for ((j=1;j<=$num; j++))
                do
                        jscript=$lmsgtf/$j.sh
                        echo "$gffcompare -M -N -r $ref -o $lmsgtf/individual.${j}.stats $lmsgtf/$j.gtf" > $jscript
                        #echo "$gtfcuff roc $cur/merge_$num/linkmerge_${num}_sgtf/individual.${j}.$j.gtf.tmap 225036 cov > $cur/merge_$num/linkmerge_${num}_sgtf/individual.${j}.roc.txt" >> $script
                        echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $lmsgtf/$j.gtf > $lmsgtf/$j.otherchrm.gtf" >> $jscript
                        echo "$gffcompare -M -N -r $ref -o $lmsgtf/$j.other.stats $lmsgtf/$j.otherchrm.gtf" >> $jscript
                        echo "$jscript" >> $lmsgtf/gff-scripts   
                        chmod u+x $jscript
                done
                
                echo "cd $lmsgtf" >> $script
                echo "cat gff-scripts | xargs -L 1 -I CMD -P 10 bash -c CMD 1> /dev/null 2> /dev/null &" >> $script
                chmod +x $script
                echo $script >> $joblist

        done

        cat $joblist | xargs -L 1 -I CMD -P 10 bash -c CMD 1> /dev/null 2> /dev/null &

        echo "Finish Aletsch - Link."
fi


# Scallop2
if [ "A" == "B" ];then
	echo "generating scallop2's assemblies..."

	sgtf=$result/scallop2
	mkdir -p $sgtf
	cd $sgtf

	joblist=$sgtf/jobs.list
	rm -rf $joblist
	for((i=1;i<=369;i++));
	do
		script=$sgtf/$i.sh
                echo "cd $sgtf" > $script 
                echo "$scallop2 -i $index/$i.bam -o $sgtf/$i.gtf > $sgtf/$i.log" >> $script
                echo "$gffcompare -M -N -r $ref -o $i.mul.stats $i.gtf" >> $script
                echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $i.gtf > $i.otherchrm.gtf" >> $script
                echo "$gffcompare -M -N -r $ref -o $i.other.stats $i.otherchrm.gtf" >> $script
		chmod +x $script
                echo $script >> $joblist
	done

	cat $joblist | xargs -L 1 -I CMD -P 20 bash -c CMD 1> /dev/null 2> /dev/null &

fi

# Stringtie2
if [ "A" == "B" ];then
        echo "generating stringtie2's assemblies..."

        sgtf=$result/stringtie2
        mkdir -p $sgtf
        cd $sgtf

        joblist=$sgtf/jobs.list
        rm -rf $joblist
        for((i=1;i<=369;i++));
        do
                script=$sgtf/$i.sh
                echo "cd $sgtf" > $script 
                echo "$stringtie2 $index/$i.bam -o $sgtf/$i.gtf > $sgtf/$i.log" >> $script
                echo "$gffcompare -M -N -r $ref -o $i.mul.stats $i.gtf" >> $script
                echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $i.gtf > $i.otherchrm.gtf" >> $script
                echo "$gffcompare -M -N -r $ref -o $i.other.stats $i.otherchrm.gtf" >> $script
                chmod +x $script
                echo $script >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 20 bash -c CMD 1> /dev/null 2> /dev/null &
fi

mkdir -p $result/cp

# copy transmeta results
if [ "A" == "B" ];then
        cur=$result/transmeta
        mkdir -p $result/cp/transmeta

        for num in $testlist;
        do
                # from gtf -> tlist, chrm and tid
                mkdir -p $result/cp/transmeta/merge_$num

                for ((j=1;j<=$num; j++));
                do
                        awk '{print $1, $12}' $cur/merge_$num/TransMeta.bam$j.gtf > $result/cp/transmeta/merge_$num/$j.tlist
                        sed -i 's/[";]//g' $result/cp/transmeta/merge_$num/$j.tlist
                        awk '{print $2, $3, $5}' $cur/merge_$num/individual.$j.TransMeta.bam$j.gtf.tmap > $result/cp/transmeta/merge_$num/$j.tmap
                done
        done
fi

# copy psiclass results
if [ "A" == "B" ];then
        cur=$result/psiclass

        rm -rf $result/cp/psiclass
        mkdir -p $result/cp/psiclass

        for num in $testlist;
        do
                # from gtf -> tlist, chrm and tid
                mkdir -p $result/cp/psiclass/merge_$num

                for ((j=1;j<=$num; j++));
                do
                        i=$((j - 1))
                        awk '{print $1, $12}' $cur/merge_$num/cov_1/psiclass_sample_$i.gtf > $result/cp/psiclass/merge_$num/$j.tlist
                        sed -i 's/[";]//g' $result/cp/psiclass/merge_$num/$j.tlist
                        awk '{print $2, $3, $5}' $cur/merge_$num/cov_1/individual.$i.psiclass_sample_$i.gtf.tmap > $result/cp/psiclass/merge_$num/$j.tmap
                done
        done
fi

# copy aletsch-link results
if [ "A" == "B" ];then
        cur=$result/aletsch-link

        rm -rf $result/cp/aletsch-link
        mkdir -p $result/cp/aletsch-link

        for num in $testlist;
        do
                # from gtf -> tlist, chrm and tid
                mkdir -p $result/cp/aletsch-link/merge_$num

                #awk '{print $2, $3, $5}' $cur/merge_$num/merge_$num.linkmerge_$num.gtf.tmap > $result/cp/aletsch-link/merge_$num/merge_$num.linkmerge_$num.gtf.tmap
                cp $cur/merge_$num/*.csv $result/cp/aletsch-link/merge_$num

                #for ((j=1;j<=$num;j++))
                #do
                        #awk '{print $1, $12}' $cur/merge_$num/linkmerge_${num}_sgtf/$j.gtf > $result/cp/aletsch-link/merge_$num/$j.tlist
                        #sed -i 's/[";]//g' $result/cp/aletsch-link/merge_$num/$j.tlist
                        #awk '{print $2, $3, $5}' $cur/merge_$num/linkmerge_${num}_sgtf/individual.$j.$j.gtf.tmap > $result/cp/aletsch-link/merge_$num/$j.tmap

                #done
        done
fi

# copy Aletsch results
if [ "A" == "B" ];then
        cur=$result/aletsch

        rm -rf $result/cp/aletsch
        mkdir -p $result/cp/aletsch

        for num in $testlist;
        do
                # from gtf -> tlist, chrm and tid
                mkdir -p $result/cp/aletsch/merge_$num

                for ((j=1;j<=$num; j++));
                do
                        awk '{print $1, $12}' $cur/merge_$num/gtf/$(($j - 1)).gtf > $result/cp/aletsch/merge_$num/$j.tlist
                        sed -i 's/[";]//g' $result/cp/aletsch/merge_$num/$j.tlist

			#$gffcompare -r $refgtf -o $cur/merge_$num/$j $cur/merge_$num/$j.gtf
                        awk '{print $2, $3, $5}' $cur/merge_$num/gtf/$(($j - 1)).mul.$(($j - 1)).gtf.tmap > $result/cp/aletsch/merge_$num/$j.tmap
                done
        done
fi
