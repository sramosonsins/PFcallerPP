#COMPILE:

#Linux
gcc *.c -lm -Wall -o PFcaller -DnoMacOS1 -lz
#MacOS
gcc *.c -lm -Wall -o PFcaller -DnoMacOS0 -lz

#RUN:

./PFcaller -h
