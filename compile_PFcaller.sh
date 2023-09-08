
#if [ "$(uname)" == "Darwin" ]; then
#function _wget() { curl "${1}" -o $(basename "${1}") ; };
#alias wget='_wget'
#fi

##zlib 1.2.11 installation (dependency)
#mkdir -p ./zlib
#wget http://zlib.net/zlib-1.2.11.tar.gz -P ./zlib
#tar -zxvf ./zlib/zlib-1.2.11.tar.gz -C ./zlib
#rm ./zlib/zlib-1.2.11.tar.gz
#cd ./zlib/zlib-1.2.11
#./configure
#make
#sudo make install

##gsl installation (dependency)
#mkdir -p /tmp/gsl
#curl -o /tmp/gsl-2.2.tar.gz ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.tar.gz -LOk
#tar -zxvf /tmp/gsl-2.2.tar.gz -C /tmp/gsl && \
#rm /tmp/gsl-2.2.tar.gz && \
#cd /tmp/gsl/gsl-2.2 && \
#./configure && \
#make && \
#sudo make install

#To compile:
gcc ./sources/*.c -lm -o ./bin/PFcaller -Wall -O3 -lz 

