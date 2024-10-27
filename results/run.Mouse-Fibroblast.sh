#!/bin/bash
dir=`pwd`
data=$dir/../data/Mouse-Fibroblast/zUMIs_output/demultiplexed
index=$dir/../data/Mouse-Fibroblast_index
ref=/data/shared/data/mouse/Mus_musculus.GRCm39.110.gtf
scallop2=$dir/../programs/scallop2
stringtie2=$dir/../programs/stringtie2
lmversion=8
linkmerge=$dir/../programs/scAletsch-${lmversion}
transmeta=$dir/../programs/TransMeta/TransMeta
psiclass=$dir/../programs/psiclass/psiclass
gffcompare=$dir/../programs/gffcompare
gtfcuff=$dir/../programs/gtfcuff
result=$dir/Mouse-Fibroblast_results

#covlist='0 0.5 1 2 5 10 20 30 50 80 100 150'
covlist='1'
testlist='2 5 10 30 50 80 100 150 200 300 369'
#testlist='10 50 100 200 300 369'
#testlist='369'

# step 0: check to-be-used tools
if [ "A" == "b" ];then
	echo "================================================================="
	echo "Check if to-be-used tools are properly linked..."
	echo "================================================================="
	if [ -e $scallop2 ];then
		echo -e "Tool Scallop2 found successfully!"
	else
		echo -e "Tool Scallop2 not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
		echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
    		exit 1
	fi

	if [ -e $stringtie2 ];then
                echo -e "Tool StringTie2 found successfully!"
        else
                echo -e "Tool StringTie2 not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi

        if [ -e $linkmerge ];then
                echo -e "Tool Linkmerge found successfully!"
        else
                echo -e "Tool Linkmerge not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi

	if [ -e $transmeta ];then
                echo -e "Tool TransMeta found successfully!"
        else
                echo -e "Tool TransMeta not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi

	if [ -e $psiclass ];then
                echo -e "Tool PsiCLASS found successfully!"
        else
                echo -e "Tool PsiCLASS not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi

	if [ -e $gffcompare ];then
                echo -e "Tool gffcompare found successfully!"
        else
                echo -e "Tool gffcompare not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi

	if [ -e $gtfcuff ];then
                echo -e "Tool gtfcuff found successfully!"
        else
                echo -e "Tool gtfcuff not found in directory 'programs'.\nPlease follow the instructions in 'Step 1: Download and Link Tools' to properly download and link all necessary tools to the directory 'programs'."
                echo -e "\nNote: Tools are not downloaded automatically. Users need to download and/or compile all required tools, and then link them to 'programs' directory before running experiments.\n"
                exit 1
        fi

	echo -e "All To-be-used tools found successfully!"

fi


#============================================
# create index for cells
#=============================================
if [ "A" == "B" ];then
        mkdir -p $index
        cd $index

        i=1
        for k in `ls $data/*.demx.bam`
        do
                if [ "$i" -lt "500"  ]; then
                        ln -sf $data/$k .
                        echo "cell #$i: $k" >> index.list
                        mv $k $i.bam
                        let i++
                fi
        done
fi


#============================================
# assembly and evaluate
#============================================
# TransMeta
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
		echo "$transmeta -B $list -s unstranded -o $cur/merge_$num > $cur/merge_$num/merge_$num.log" > $script
                echo "cp $cur/merge_$num/TransMeta.gtf $cur/merge_$num/TransMeta-2.gtf" >> $script

		#for cov in $covlist;
		#do
			#echo "$gffcompare -r $ref -o $cur/merge_$num/merge_${num}_$cov.stats $cur/merge_$num/TransMeta-$cov.gtf" >> $script
			#echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $cur/merge_$num/TransMeta-$cov.gtf > $cur/merge_$num/TransMeta-$cov-otherchrm.gtf" >> $script
                        #echo "$gffcompare -r $ref -o $cur/merge_$num/merge_${num}_$cov-otherchrm.stats $cur/merge_$num/TransMeta-$cov-otherchrm.gtf" >> $script
		#done

		# gffcompare individual-level;1-index in TransMeta
                for ((j=1;j<=$num; j++))
                do
                        echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/individual.${j}.stats $cur/merge_$num/TransMeta.bam$j.gtf" >> $script
                        echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $cur/merge_$num/TransMeta.bam$j.gtf > $cur/merge_$num/TransMeta.bam$j.otherchrm.gtf" >> $script
                        echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/individual.${j}.other.stats $cur/merge_$num/TransMeta.bam$j.otherchrm.gtf" >> $script
                done

		chmod +x $script
		echo $script >> $joblist
	done

	cat $joblist | xargs -L 1 -I CMD -P 12 bash -c CMD 1> /dev/null 2> /dev/null &

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

		for cov in $covlist;
		do
			mkdir -p $cur/merge_${num}/cov_$cov
                	echo "$psiclass --lb $list -o $cur/merge_$num/cov_$cov/ --vd $cov -p 10 > $cur/merge_$num/cov_$cov/merge_${num}_$cov.log" >> $script
                        #echo "$gffcompare -r $ref -o $cur/merge_$num/merge_${num}_$cov.stats $cur/merge_$num/cov_$cov/psiclass_vote.gtf" >> $script
                	#echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $cur/merge_$num/cov_$cov/psiclass_vote.gtf > $cur/merge_$num/cov_$cov/psiclass_vote-otherchrm.gtf" >> $script
                        #echo "$gffcompare -r $ref -o $cur/merge_$num/merge_${num}_$cov-otherchrm.stats $cur/merge_$num/cov_$cov/psiclass_vote-otherchrm.gtf" >> $script
		done

		# gffcompare individual-level;0-index in PsiClass
                for ((j=0;j<$num; j++))
                do
                        echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/cov_1/individual.${j}.stats $cur/merge_$num/cov_1/psiclass_sample_$j.gtf" >> $script
                        echo "awk '\$1 !~ /^(1|2|3|4|5|6|7|8|9)$/' $cur/merge_$num/cov_1/psiclass_sample_$j.gtf > $cur/merge_$num/cov_1/psiclass_sample_$j.otherchrm.gtf" >> $script
                        echo "$gffcompare -M -N -r $ref -o $cur/merge_$num/individual.${j}.other.stats $cur/merge_$num/cov_1/psiclass_sample_$j.otherchrm.gtf" >> $script
                done

                chmod +x $script
                echo $script >> $joblist
        done

        cat $joblist | xargs -L 1 -I CMD -P 12 bash -c CMD 1> /dev/null 2> /dev/null &

        echo "Finish PsiCLASS."
fi

# Aletsch
aletschExe=aletsch-585c8a 
aletsch=$dir/../programs/$aletschExe
model=$dir/../models/aletsch-39.ensembl.chr1.joblib

if [ "A" == "B" ];then
        echo "Running Aletsch..."

        cur=$result/aletsch
        #rm -rf $cur
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
                    echo "$index/$i.bam $index/$i.bai paired_end" >> $list
                done

                test=$cur/merge_$num
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
                        echo "$gffcompare -M -N -r $ref -o $k.other.stats $k.otherchrm.gtf" >> $kscript
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

# copy aletsch-link results
if [ "A" == "B" ];then
        cur=$result/aletsch-link

        rm -rf $result/cp/aletsch-link
        mkdir -p $result/cp/aletsch-link

        for num in $testlist;
        do
                # from gtf -> tlist, chrm and tid
                mkdir -p $result/cp/aletsch-link/merge_$num

                awk '{print $2, $3, $5}' $cur/merge_$num/merge_$num.linkmerge_$num.gtf.tmap > $result/cp/aletsch-link/merge_$num/merge_$num.linkmerge_$num.gtf.tmap
                cp $cur/merge_$num/*.csv $result/cp/aletsch-link/merge_$num

                for ((j=1;j<=$num;j++))
                do
                        awk '{print $1, $12}' $cur/merge_$num/linkmerge_${num}_sgtf/$j.gtf > $result/cp/aletsch-link/merge_$num/$j.tlist
                        sed -i 's/[";]//g' $result/cp/aletsch-link/merge_$num/$j.tlist
                        awk '{print $2, $3, $5}' $cur/merge_$num/linkmerge_${num}_sgtf/individual.$j.$j.gtf.tmap > $result/cp/aletsch-link/merge_$num/$j.tmap

                done
        done
fi