################ Exploratory Data Analysis Part 2 - simplified version #############################


# In part 2 of EDA, we are looking at mapping statistics and how the reads compare before and after removing
# PCR duplicates

#setup

library(here) #this should set your working directory to your project folder where this script is located
list.files()

library(ggplot2)
library(tidyr)
library(data.table)
library(reshape2)
#Multiplot Functions####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##### Read in files and reformat

list.files(path="../",pattern="All_head") #list files with this pattern in the title
mapped.raw<-read.delim("../All_head_flagstat_reformat.txt",header=T) #read in file generated after mapping
mapped.fltr<-read.delim("../All_head_filt_flagstat_reformat.txt",header=T) #read in file generated after removing PCR dups
mapped.raw$ID<-gsub("_sort_flagstat.txt", "", mapped.raw$sampleID)
mapped.fltr$ID<-gsub("_sortfltr_flagstat.txt", "", mapped.fltr$sampleID)
mapped.comb<-merge(mapped.raw,mapped.fltr, by="ID") #make one dataframe with both mapped read types
mapped.comb<-mapped.comb[,c(1,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22,11,23,12,24,13,25)] #reorganize dataframe, this can be changed if needed, but this column format generally works
vars<-colsplit(mapped.comb$ID, "_", c("plate","LABID","mpatch","year","well")) #the options after c() will be specifiy to your naming scheme for your samples, CHANGE FOR YOUR DATA
mapped.comb<-cbind(mapped.comb,vars)
mapped.comb<-mapped.comb[,c(1,24:28,2:23)] #tidy up


#get proportions of filtered to unfiltered mapped reads
mapped.comb$prop.NDreads<-mapped.comb$mapped.reads.QC.passed.y/mapped.comb$total.reads.QC.passed.x
summary(mapped.comb$prop.NDreads)#see percentage of raw reads left after filtering
sd(mapped.comb$prop.NDreads) #get sd of how many reads are dropped to filtering
hist(mapped.comb$prop.NDreads) #plot histogram of prop of nonduplicated reads
mapped.comb$prop.dupreads<-(1-mapped.comb$prop.NDreads) #get proportion of reads that were duplicates
summary(mapped.comb$prop.dupreads)
sd(mapped.comb$prop.dupreads)
hist(mapped.comb$prop.dupreads)

mapped.comb.short<-mapped.comb[,c(1:12,29)] #create shorter version of he above dataframe
str(mapped.comb.short)
names(mapped.comb.short) <- c("ID","plate","LABID","mpatch","year","well", "UF_total_reads_QC_passed"	,"F_total_reads_QC_passed","UF_mapped_reads_QC_passed","F_mapped_reads_QC_passed",	"UF_percent_reads_mapped","F_percent_reads_mapped","prop") #you should change some of these if the variables aren't in your dataset, like mpatch
hist(mapped.comb.short$F_percent_reads_mapped)
write.csv(mapped.comb.short,"../short_mapping_stats.csv") #change path to where you want your mapping stat file written to

#add in QC category info from section I:
# you can do this step if you included your failed sampels during your mapping step, if not this can be skipped. Basically this just allows us to look at EDA with only qc-passed samples
SRA1<-sample.reads.binRAshort[,c(1,5:9)] #this dataframe was generated in EDA part 1 script, which is upstream of this one. Can copy from that to get qc info.
SRA1$platewell<-paste(SRA1$plate,SRA1$well)
mapped.comb.short$platewell<-paste(mapped.comb.short$plate,mapped.comb.short$well)
CFS<-merge(mapped.comb.short,SRA1,by="platewell")
CFS.ok<-subset(CFS,QC_cat!="failed")

#Input concentration matters, check effects of input concentration
# To do this, you need the total input DNA for each well, either in its own sheet like below or as a column in your metadata sheet to pull in

DNAconc<-read.csv("inputDNAconckey.csv") #this should specify the file that has the input DNA concentrations for your library prep per sample
CFS.ok1<-merge(DNAconc,CFS.ok)#if you skpped above step, can just merge DNAconc with the mapped.comb.short dataframe
CFS.ok1$concentration<-factor(CFS.ok1$concentration, levels=c("<=50","<=75 & >=50","100"))# can also try skipping this and just plotting with all concentrations
DNAbox<-ggplot(CFS.ok1, aes(x=concentration, y=prop)) +theme_bw()+
  geom_boxplot(colour="black", fill="red", alpha=0.7)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4)+
  ylab("Proportion Filtered Mapped/Sequenced Fragments") +xlab("DNA Input Original Concentration (ng/ul)")

pdf("prop filtered mapped-seq frag by DNA input.pdf", width=6, height=6)
dev.off()


