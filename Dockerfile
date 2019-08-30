FROM ubuntu:16.04

RUN apt-get update
RUN apt install -y zlib1g-dev gfortran g++ make curl liblzma-dev libpcre3 libpcre3-dev libcurl4-openssl-dev

RUN apt-get update
RUN apt install -y libssl-dev libxml2-dev software-properties-common sysstat libbz2-dev vim wget git

RUN add-apt-repository ppa:ubuntugis/ppa
RUN apt-get update
RUN apt-get install -y gdal-bin libgdal1-dev

WORKDIR /home/saps
COPY . .

#RUN ./R-3.6.1/configure --with-readline=no --with-x=no
#RUN cd R-3.6.1 && make && make install
#RUN make
#RUN make install R-3.6.1
