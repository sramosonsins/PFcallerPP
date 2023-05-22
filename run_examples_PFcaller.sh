#To compile:
gcc ./sources/*.c -lm -o ./bin/PFcaller -Wall -O3 -lz
cd ./Examples

echo
echo flags included in Pfcaller_help.txt
../bin/PFcaller -S ./chr1 -h > ../Pfcaller_help.txt 

#Run examples
echo --------------------------------------------------------------------------------------------------
echo Toy file: 10bp. TESTING OPTIONAL PARAMETERS WITH ploidy 1
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex10bpp01.01.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f t -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f t -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 

FILE_RES="../Results/pileup_test00p01.tfa"
FILE_VAL="../Validation/pileup_test00p01.tfa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo Example Ex10bpp01.0b.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f f -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f f -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 

FILE_RES="../Results/pileup_test00p01_chr1.fa"
FILE_VAL="../Validation/pileup_test00p01_chr1.fa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo Example Ex10bpp01.01c.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 

FILE_RES="../Results/pileup_test00p01_lik.txt"
FILE_VAL="../Validation/pileup_test00p01_lik.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/pileup_test00p01_complete_table.txt"
FILE_VAL="../Validation/pileup_test00p01_complete_table.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo Example Ex10bpp01.01d.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 

FILE_RES="../Results/pileup_test00p01.gvf"
FILE_VAL="../Validation/pileup_test00p01.gvf"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo --------------------------------------------------------------------------------------------------
echo Toy file: 10bp. TESTING PARAMETERS WITH ploidy 2
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex10bpp02.01.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 2  -s 13743 -m 5 -M 20 -q 20 -r 1 -f 0vgft -i ./ex16.pileup -o ../Results/pileup_test00p02
../bin/PFcaller -S ./chr1_20.txt -p 2  -s 13743 -m 5 -M 20 -q 20 -r 1 -f 0vgft -i ./ex16.pileup -o ../Results/pileup_test00p02 

FILE_RES="../Results/pileup_test00p02.tfa"
FILE_VAL="../Validation/pileup_test00p02.tfa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/pileup_test00p02.gvf"
FILE_VAL="../Validation/pileup_test00p02.gvf"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/pileup_test00p02_chr1.fa"
FILE_VAL="../Validation/pileup_test00p02_chr1.fa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/pileup_test00p02_lik.txt"
FILE_VAL="../Validation/pileup_test00p02_lik.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/pileup_test00p02_complete_table.txt"
FILE_VAL="../Validation/pileup_test00p02_complete_table.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo --------------------------------------------------------------------------------------------------
echo Toy file: 10bp. TESTING OPTIONAL PARAMETERS WITH ploidy 10
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex10bpp10.01.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f t -i ./ex16b.pileup -o ../Results/pileup_test00p10
../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f t -i ./ex16b.pileup -o ../Results/pileup_test00p10 

FILE_RES="../Results/pileup_test00p10.tfa"
FILE_VAL="../Validation/pileup_test00p10.tfa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES


echo 
echo Example Ex10bpp10.01b.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f f -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10
../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f f -i ./ex16b.pileup -o ../Results/pileup_test00p10 

FILE_RES="../Results/pileup_test00p10_chr1.fa"
FILE_VAL="../Validation/pileup_test00p10_chr1.fa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo Example Ex10bpp10.01c.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10
../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f 0 -i ./ex16b.pileup -o ../Results/pileup_test00p10 

FILE_RES="../Results/pileup_test00p10_lik.txt"
FILE_VAL="../Validation/pileup_test00p10_lik.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/pileup_test00p10_complete_table.txt"
FILE_VAL="../Validation/pileup_test00p10_complete_table.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo Example Ex10bpp10.01d.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10
../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f g -i ./ex16b.pileup -o ../Results/pileup_test00p10 

FILE_RES="../Results/pileup_test00p10.gvf"
FILE_VAL="../Validation/pileup_test00p10.gvf"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo 
echo --------------------------------------------------------------------------------------------------
echo Toy file: 20bp. TESTING OPTIONAL PARAMETERS WITH ploidy 16
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex16bpp16.01.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -i ./ex16.pileup -o ../Results/ex16
../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -i ./ex16.pileup -o ../Results/ex16

FILE_RES="../Results/ex16.tfa"
FILE_VAL="../Validation/ex16.tfa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo
#echo Example Ex16bpp16.02.txt 
#echo ../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex16b
#../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex16b 
#echo
echo --------------------------------------------------------------------------------------------------
echo Toy file: 20bp. TESTING OPTIONAL PARAMETERS WITH ploidy 1000000
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex16bpp1e6.01.txt ---------------------------------
#echo ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 20 -q 20 -f t0gf -r 1 -i ./ex160p.pileup -o ../Results/ex1e6p
../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 20 -q 20 -f t0gf -r 1 -i ./ex16b.pileup -o ../Results/ex1e6p

FILE_RES="../Results/ex1e6p.tfa"
FILE_VAL="../Validation/ex1e6p.tfa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p.gvf"
FILE_VAL="../Validation/ex1e6p.gvf"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p_complete_table.txt"
FILE_VAL="../Validation/ex1e6p_complete_table.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p_lik.txt"
FILE_VAL="../Validation/ex1e6p_lik.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p_lik.txt"
FILE_VAL="../Validation/ex1e6p_lik.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p_lik.txt"
FILE_VAL="../Validation/ex1e6p_lik.txt"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES


FILE_RES="../Results/ex1e6p_chr1.fa"
FILE_VAL="../Validation/ex1e6p_chr1.fa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p_chr2.fa"
FILE_VAL="../Validation/ex1e6p_chr2.fa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

FILE_RES="../Results/ex1e6p_chr3.fa"
FILE_VAL="../Validation/ex1e6p_chr3.fa"

I=`wc -c $FILE_RES\.gz | awk '{print $1}'`
J=`wc -c $FILE_VAL\.gz | awk '{print $1}'`
if [ $I -ne $J ]; then 
  echo "****** WARNING! THE FILE $FILE_RES.gz HAS DIFFERENT SIZE THAN VALIDATION FILE! ******" ; 
fi

gunzip -k $FILE_RES\.gz
K=`diff -q $FILE_RES $FILE_VAL | awk '{print $1}'`
if [ $K ]; then
  echo "****** WARNING! THE FILE $FILE_RES HAS DIFFERENT CONTENT THAN VALIDATION FILE! ******";
  diff -ia $FILE_RES $FILE_VAL;
fi
rm $FILE_RES

echo
#echo Example Ex16bpp1e6.02.txt
#echo ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb
#../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb 
#echo
#echo Example Ex16bpp1e6.03.txt
#echo ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 5 -q 20 -f g -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb
#../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 5 -M 5 -q 20 -f g -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb 
#echo
exit;



echo --------------------------------------------------------------------------------------------------
echo 100Kb file: 200bp. TESTING OPTIONAL PARAMETERS WITH ploidy 10000
echo --------------------------------------------------------------------------------------------------
echo
echo Example LargeReadDepthFile.txt
#echo ../bin/PFcaller -S ./chr1_20.txt0 -p 10000 -s 12345 -m 10 -M 4000 -q 20 -f t -r 1 -i ./test.iSNPcall.replicate1.pool.pileup -o ../Results/test.iSNPcall.replicate1.pool
../bin/PFcaller -S ./chr1_20.txt0 -p 10000 -s 12345 -m 10 -M 4000 -q 20 -f f -r 1 -i ./test.iSNPcall.replicate1.pool.pileup -o ../Results/test.iSNPcall.replicate1.pool 
echo
echo --------------------------------------------------------------------------------------------------
echo 100Kb file: TESTING OPTIONAL PARAMETERS WITH ploidy 2
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex100Kbp02.01.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 47743 -m 5 -M 20 -q 20 -f f -i ./test.iSNPcall.replicate2.Ind1.pileup -o ../Results/test.iSNPcall.replicate2.Ind1.pileup_p02
../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 47743 -m 5 -M 20 -q 20 -f f -i ./test.iSNPcall.replicate2.Ind1.pileup -o ../Results/test.iSNPcall.replicate2.Ind1.pileup_p02 
echo 
echo Example Ex100Kbp02.03.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 57568 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate4.Ind6.pileup -o ../Results/test.iSNPcall.replicate4.Ind6.pileup_p02
../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 57568 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate4.Ind6.pileup -o ../Results/test.iSNPcall.replicate4.Ind6.pileup_p02 
echo 
echo Example Ex100Kbp02.04.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 66692 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1B.Ind1.pileup -o ../Results/test.iSNPcall.replicate1B.Ind1.pileup_p02
../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 66692 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1B.Ind1.pileup -o ../Results/test.iSNPcall.replicate1B.Ind1.pileup_p02 
echo 
#echo Example Ex100Kbp02.05.txt
echo ../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 33794 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p02
../bin/PFcaller -S ./chr1_100000.txt -p 2  -s 33794 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p02 
echo 
echo --------------------------------------------------------------------------------------------------
echo 100Kb file: TESTING OPTIONAL PARAMETERS WITH ploidy 10
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex100Kbp10.01.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 76665 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate2.Ind1.pileup -o ../Results/test.iSNPcall.replicate2.Ind1.pileup_p10
../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 76665 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate2.Ind1.pileup -o ../Results/test.iSNPcall.replicate2.Ind1.pileup_p10 
echo 
echo Example Ex100Kbp10.03.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 43886 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate4.Ind6.pileup -o ../Results/test.iSNPcall.replicate4.Ind6.pileup_p10
../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 43886 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate4.Ind6.pileup -o ../Results/test.iSNPcall.replicate4.Ind6.pileup_p10 
echo 
echo Example Ex100Kbp10.04.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 86231 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1B.Ind1.pileup -o ../Results/test.iSNPcall.replicate1B.Ind1.pileup_p10
../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 86231 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1B.Ind1.pileup -o ../Results/test.iSNPcall.replicate1B.Ind1.pileup_p10 
echo 
echo Example Ex100Kbp10.05.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 68923 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p10
../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 68923 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p10 
echo 
echo Example Ex100Kbp10.06.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 68923 -m 5 -M 20 -q 20 -f 0 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p10
../bin/PFcaller -S ./chr1_100000.txt -p 10 -s 68923 -m 5 -M 20 -q 20 -f 0 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p10 
echo

echo --------------------------------------------------------------------------------------------------
echo 100Kb file: TESTING OPTIONAL PARAMETERS WITH ploidy 1000000
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex100Kbp1e6.01.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 76665 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate2.Ind1.pileup -o ../Results/test.iSNPcall.replicate2.Ind1.pileup_p1e6
../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 76665 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate2.Ind1.pileup -o ../Results/test.iSNPcall.replicate2.Ind1.pileup_p1e6 
echo
echo Example Ex100Kbp1e6.03.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 43886 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate4.Ind6.pileup -o ../Results/test.iSNPcall.replicate4.Ind6.pileup_p1e6
../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 43886 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate4.Ind6.pileup -o ../Results/test.iSNPcall.replicate4.Ind6.pileup_p1e 
echo 
echo Example Ex100Kbp1e6.04.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 86231 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1B.Ind1.pileup -o ../Results/test.iSNPcall.replicate1B.Ind1.pileup_p1e6
../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 86231 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1B.Ind1.pileup -o ../Results/test.iSNPcall.replicate1B.Ind1.pileup_p1e6 
echo 
echo Example Ex100Kbp1e6.05.txt
#echo ../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 68923 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p1e6
../bin/PFcaller -S ./chr1_100000.txt -p 1000000 -s 68923 -m 5 -M 20 -q 20 -i ./test.iSNPcall.replicate1.Ind1.pileup -o ../Results/test.iSNPcall.replicate1.Ind1.pileup_p1e6 
echo 
