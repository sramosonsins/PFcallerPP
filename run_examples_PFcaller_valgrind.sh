#To compile:
gcc ./sources/*.c -lm -o ./bin/PFcaller -Wall -lz -g
cd ./Examples

echo
echo flags included in Pfcaller_help.txt
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S chr1 -h > ../Pfcaller_help.txt 2>../Results/valgrind00.txt

#Run examples
echo --------------------------------------------------------------------------------------------------
echo Toy file: 10bp. TESTING OPTIONAL PARAMETERS WITH ploidy 1
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex10bpp01.01.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f t -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f t -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 2>../Results/valgrind01.txt
echo 
echo Example Ex10bpp01.0b.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f f -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f f -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 2>../Results/valgrind01b.txt
echo 
echo Example Ex10bpp01.01c.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 2>../Results/valgrind01c.txt
echo 
echo Example Ex10bpp01.01d.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 1  -s 12345 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p01 2>../Results/valgrind01d.txt
echo 
echo --------------------------------------------------------------------------------------------------
echo Toy file: 10bp. TESTING PARAMETERS WITH ploidy 2
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex10bpp02.01.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 2  -s 13743 -m 5 -M 20 -q 20 -r 1 -i ./ex16.pileup -o ../Results/pileup_test00p02
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 2  -s 13743 -m 5 -M 20 -q 20 -r 1 -i ./ex16.pileup -o ../Results/pileup_test00p02 2>../Results/valgrind02.txt
echo 
echo Example Ex10bpp02.01.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 2  -s 13743 -m 5 -M 20 -q 20 -r 1 -i ./ex16.pileup.gz -o ../Results/pileup_test00p02
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 2  -s 13743 -m 5 -M 20 -q 20 -r 1 -i ./ex16.pileup.gz -o ../Results/pileup_test00p02 2>../Results/valgrind02gz.txt
echo 
echo --------------------------------------------------------------------------------------------------
echo Toy file: 10bp. TESTING OPTIONAL PARAMETERS WITH ploidy 10
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex10bpp10.01.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f t -i ./ex16b.pileup -o ../Results/pileup_test00p10
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f t -i ./ex16b.pileup -o ../Results/pileup_test00p10 2>../Results/valgrind04.txt
echo 
echo Example Ex10bpp10.01b.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f f -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -r 1 -f f -i ./ex16b.pileup -o ../Results/pileup_test00p10 2>../Results/valgrind04b.txt
echo 
echo Example Ex10bpp10.01c.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f 0 -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10 2>../Results/valgrind04c.txt
echo 
echo Example Ex10bpp10.01d.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr1_20.txt -p 10 -s 75653 -m 5 -M 20 -q 20 -f g -r 1 -i ./ex16b.pileup -o ../Results/pileup_test00p10 2>../Results/valgrind04d.txt
echo 
echo --------------------------------------------------------------------------------------------------
echo Toy file: 20bp. TESTING OPTIONAL PARAMETERS WITH ploidy 16
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex16bpp16.01.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -i ./ex16.pileup -o ../Results/ex16
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -i ./ex16.pileup -o ../Results/ex16 2>../Results/valgrind06.txt
echo
#echo Example Ex16bpp16.02.txt 
#echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex16b
#valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 16  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex16b 2>../Results/valgrind07.txt
#echo
echo --------------------------------------------------------------------------------------------------
echo Toy file: 20bp. TESTING OPTIONAL PARAMETERS WITH ploidy 1000000
echo --------------------------------------------------------------------------------------------------
echo 
echo Example Ex16bpp1e6.01.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -i ./ex16.pileup -o ../Results/ex1e6p
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -i ./ex16.pileup -o ../Results/ex1e6p 2>../Results/valgrind08.txt
echo
echo Example Ex16bpp1e6.02.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 5 -M 5 -q 20 -f t -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb 2>../Results/valgrind09.txt
echo
echo Example Ex16bpp1e6.02.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 5 -q 20 -f g -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 5 -M 5 -q 20 -f g -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6bb 2>../Results/valgrind09b.txt
echo
echo Example Ex16bpp1e6.02.txt
echo valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 2 -M 400 -q 20 -f g -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6c
valgrind --leak-check=full --track-origins=yes -v ../bin/PFcaller -S ./chr12_2030.txt -p 1000000  -s 12345 -m 5 -M 400 -q 20 -f g -r 1 -B 0 -i ./ex16b.pileup -o ../Results/ex1e6c 2>../Results/valgrind09c.txt
echo
