library('getopt');

params=c(
	"input_file", "i", 1, "character",
	"output_file", "o", 2, "character"
);

opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE), debug=FALSE);
script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2];


InputFileName=opt$input_file;
OutputRoot=gsub(".summary_table.xls", "", InputFileName);

if(length(opt$output_file)){
	OutputRoot=opt$output_file;
}

cat("Input File: ", InputFileName, "\n");
cat("Output File Root: ", OutputRoot, "\n");

pdf(paste(OutputRoot, ".pdf", sep=""), height=11, width=8.5);

################################################################################

load_summary_table=function(filename){
	st=as.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1, check.names=F));
	return(st[,2:ncol(st)]);
}

normalize=function(st){
	sums=apply(st, 1, sum);
	n=matrix(0,nrow=nrow(st), ncol=ncol(st));
	rownames(n)=rownames(st);
	colnames(n)=colnames(st);
	for(i in 1:nrow(st)){
		n[i,]=st[i,]/sums[i];
	}
	return(n);
}

################################################################################

# Load the summary table
st=load_summary_table(InputFileName);
num_taxa=ncol(st);
num_samples=nrow(st);
taxa_names=colnames(st);
sample_names=rownames(st);

cat("Num Taxa: ", num_taxa, "\n");
cat("Num Samples: ", num_samples, "\n");
#print(st);

# Normalize counts
nst=normalize(st);
#print(nst);

mean_prop=numeric(num_taxa);
sd_prop=numeric(num_taxa);
valid=logical(num_taxa);
ubiquity=numeric(num_taxa);

# For each taxa, compute the ubiquity and mean(log(proportions))
for(i in 1:num_taxa){
	prop=nst[,i];
	prop=prop[log(prop,10)>-4];	# Treat as non response
	n=length(prop);
	taxa=taxa_names[i];

	ubiquity[i]=n/num_samples;

	prop=log(prop,10);

	cat("Taxa: ", taxa, "\n");

	print(prop);
	mean_prop[i]=mean(prop);
	sd_prop[i]=sd(prop);

	if(ubiquity[i]>=.05 && n>1){
		cat("  Mean(Log[p]):", mean_prop[i], "\n");
		cat("  SD(Log[p]):", sd_prop[i], "\n");
		valid[i]=TRUE;
	}else{
		valid[i]=FALSE;
	}		
}

# Subset out the taxa that are can be drawn
mean_prop=mean_prop[valid];
sd_prop=sd_prop[valid];
taxa_names=taxa_names[valid];
mean_sd_cor=cor(mean_prop, sd_prop);
ubiquity=ubiquity[valid];
num_valid=sum(valid);

# Estimate correlation between var and mean
#r=mean_sd_cor;
#t=r*(sqrt(num_valid-2)/sqrt(1-r^2));
#cor_pval=1-pt(t,num_valid-2);

# Compute rank of outlierness
#reg_line=lm(sd_prop~mean_prop);
#outlier_rank=num_valid-rank(abs(reg_line$residuals))+1;

# Plot taxa variances
ylim=range(sd_prop);
xlim=range(mean_prop);

plot(mean_prop, sd_prop, 
	xlim=c(xlim[1]*1.05, .5),
	ylim=c(ylim[1]*.8, ylim[2]*1.05),
	xlab="Mean(Log[Abundance])", ylab="St. Dev.(Log[Abundance])", main="Variation vs. Abundance", type="n");

# Label plot with input file name
mtext(OutputRoot, line=0, font=2, cex=.7);

#mtext(sprintf("Correlation: %3.2f", mean_sd_cor), line=-1);
#mtext(sprintf("p-value: %.3f", cor_pval), line=-2);

# Draw regression line
#abline(reg_line, col="grey", lwd=5);

# Plot point for each taxa
points(mean_prop, sd_prop, cex=4*sqrt(ubiquity/pi));

# Label taxa
dist_mat=as.matrix(dist(cbind(mean_prop, sd_prop)), method="manhattan");
max_dist=max(dist_mat);
dist_mat=dist_mat+diag(max_dist,num_valid);
min_dist=apply(dist_mat, 1, min);

# This will label the top 10 outliers with ubiquity>.05, with mean(log[abundance])>=-3
#	label_ix=(ubiquity>.05 && outlier_rank<10) | mean_prop>=-3;

# If you only want to label certain taxa, do the following here:
#	label_ix=(taxa_names=="Bacilli_unclassified" | taxa_names=="Actinomycetales_unclassified" | taxa_names=="Peptostreptococcaceae_Peptostreptococcus")

# This will label all the taxa that are the most separated (min_dist>.07) or with significant mean abundance
label_ix=(min_dist>.10 & mean_prop>-3.2) | (min_dist>.06 & mean_prop >-2)

text(mean_prop[label_ix], sd_prop[label_ix], label=taxa_names[label_ix], col="darkred", pos=3, cex=1.2);

# Legend
cutoffs=c(.05, .25, .5, .75, 1);
labels=sprintf("%2.0f%%", 100*cutoffs);
sizes=4*sqrt(cutoffs/pi)
legend(0, ylim[2]*1.04, legend=c(labels), pch=1, pt.cex=c(sizes), title="Ubiquity");

################################################################################

cat("Done.\n")

dev.off();
q(status=0)
