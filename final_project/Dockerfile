FROM ubuntu
RUN apt-get update
RUN apt-get -y install openjdk-8-jre g++ make wget git bwa samtools

WORKDIR /tmp
RUN git clone git://github.com/GregoryFaust/samblaster.git
WORKDIR /tmp/samblaster
RUN make
RUN cp samblaster /usr/local/bin/.

WORKDIR /opt
RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
RUN tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2
ENV PATH "$PATH:/opt/strelka-2.9.10.centos6_x86_64/bin"
RUN wget https://github.com/Illumina/manta/releases/download/v1.5.0/manta-1.5.0.centos6_x86_64.tar.bz2
RUN tar xvjf manta-1.5.0.centos6_x86_64.tar.bz2
ENV PATH "$PATH:/opt/manta-1.5.0.centos6_x86_64/bin"

WORKDIR /opt
RUN mkdir /VarScan
WORKDIR /opt/VarScan
RUN wget https://nchc.dl.sourceforge.net/project/varscan/VarScan.v2.3.9.jar
RUN ln -s  VarScan.v2.3.9.jar VarScan.jar
ENV VARSCAN_PATH /opt/VarScan

RUN apt-get update \
        && apt-get -y install gnuplot \
        && DEBIAN_FRONTEND=noninteractive apt-get -y install r-base

COPY DNAcopy_install.R /tmp
RUN Rscript /tmp/DNAcopy_install.R

WORKDIR /home
