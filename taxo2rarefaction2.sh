


#!/bin/bash


## Example :
# sh taxo2rarefaction2.sh `find . -name *.CFQD_clean.greedym1_full.kaiju.full_labels.virus.blastn.best_score.taxo.virus.xls`

for i in $@ ;do sample=$(echo `basename $i` | sed -e "s/\..*$//"); echo $sample; cut -f 13 $i | sort | uniq -c | sort -k1,1nr | sed -e "s/^ *//" | awk 'BEGIN{sum=0};{sum+=$1; printf("%d %d ", $1, sum); for (ii=2;ii<=NF;ii++){printf ("%s ", $ii);}; printf("\n")};' | cut -d" " -f 2 | nl > /tmp/taxo_${sample}.tmp  ; done

y_max=$(cat /tmp/taxo_*.tmp | sort -k 2,2n | tail -n 1 | cut -f 2)
x_max=$(wc -l /tmp/taxo_*.tmp | grep tmp | sort -k 1,1n | tail -n 1 | awk '{print $1}')

~/SoftWare/R-3.3.0/R-3.3.0/BUILD/bin/R --vanilla --slave <<EEE

rm(list = ls())      # Clear all variables
graphics.off()    # Close graphics windows


mylist <- list.files("/tmp/", pattern="taxo_.*.tmp$")

mylist2 <- gsub ("taxo_", "", mylist)
mylist2 <- gsub (".tmp", "", mylist2)



plot_colors <-seq(1:length(mylist));

png("/tmp/plot.png")

par(oma=c(0,0,0,5))

x <- read.table(paste("/tmp/", mylist[1],sep=""))

plot(x\$V1, x\$V2, xlab="Total Number of Taxonomies Identified for y Sequences", ylim = c(0,$y_max), xlim=c(0,$x_max), ylab="Cumulative Number of Sequences", main="Rarefaction Curves", type="l", col=plot_colors[1]);


for (k in 2:length(mylist)){

x <- read.table(paste("/tmp/", mylist[k], sep=""))

lines(x\$V1, x\$V2, type="l", col=plot_colors[k]);

}
legend(par('usr')[2]+1, par('usr')[4]+1, xpd=NA, mylist2, col=plot_colors, lty=1,pch=NA, pt.cex=1, cex=0.5  )
dev.off()
EEE
