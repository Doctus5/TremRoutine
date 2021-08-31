#!/bin/sh

cd programs/NonLinLoc/
./clean.sh
cd ..
cd ..

jsonfiles=`basename -a results/events/*.json*`

for file in $jsonfiles
do
	echo $file
	python NNLarange.py $file 0
	cd programs/NonLinLoc/
	./NLLoc run/nlloc_LUSI.in
	cd ..
	cd ..
	python NNLarange.py $file 1
	cd programs/NonLinLoc/
	./clean.sh
	cd ..
	cd ..
done
