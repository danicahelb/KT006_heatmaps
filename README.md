# KT006_heatmaps

```r
rm(list=ls())
setwd("~/Desktop/Microarray Analyses/Kevin's array/151118 data/")
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
```

##########################################################################################
##Import clinical data
```r
print(load("./Data/Clinical Data/clinical_Apac_Jan2012.RData"))
clinical <- data
rm(data)
```

##########################################################################################
##Import array data 
```r
print(load("./Data/Processed Data/0004.KTarray.Apac.transformed_normalized.RData"))
```

##########################################################################################
##change names to follow PNAS code
```r
apac_spots <- annAnti
rm(annAnti)

pheSera$idno <- gsub("b", "", pheSera$Sample)
apac_clinical <- merge(clinical, pheSera, by="idno")
names(apac_clinical)[names(apac_clinical)=="SampleID"] <- "sampleID"
rm(pheSera, clinical)

rownames(intensity.MedA) <- gsub("15-34", "1534", rownames(intensity.MedA))
rownames(intensity.MedA) <- gsub("15-38", "1538", rownames(intensity.MedA))
rownames(intensity.MedA) <- gsub("17-34", "1734", rownames(intensity.MedA))
rownames(intensity.MedA) <- gsub("52-42", "5242", rownames(intensity.MedA))
```

####################################################################################################
#DETERMINE INTENSITY OF RESPONSE FOR EACH SAMPLE:
##for each sample, what is the mean intensity to all spots?
```r
expr <- intensity.MedA
dim(expr)
intensity.1 <- data.frame("sampleID"=colnames(expr), "mean_intensity"=apply(expr, 2, mean, na.rm=T))
head(intensity.1)
table(intensity.1$mean_intensity)
apac_clinical <- merge(apac_clinical, intensity.1, by="sampleID", all.x=T, all.y=T)
dim(apac_clinical)
#[1] 990 153
```

####################################################################################################
#DETERMINE MEAN INTENSITY OF RESPONSE FOR EACH SPOT:
```r
temp <- data.frame("mean_intensity"=apply(intensity.MedA[,colnames(intensity.MedA) %in% apac_clinical$sampleID],1,mean, na.rm=T))
temp$Unique.ID <- rownames(temp)
apac_spots <- merge(apac_spots, temp, by="Unique.ID", all.x=T)
```

##############################################################################
##create UID so Unique.ID in apac_spots to match intensity data after transposed
```r
head(apac_spots)
apac_spots$UID <- apac_spots$Unique.ID
apac_spots$UID <- gsub("^1", "X1", apac_spots$UID)
apac_spots$UID <- gsub("-", ".", apac_spots$UID)
apac_spots$UID <- gsub("^2", "X2", apac_spots$UID)
apac_spots$UID <- gsub("^3", "X3", apac_spots$UID)
apac_spots$UID <- gsub("^5", "X5", apac_spots$UID)
apac_spots$UID <- gsub("^6", "X6", apac_spots$UID)
```

##############################################################################
##save
```r
save(apac_clinical, apac_spots, intensity.MedA, intensity.MNA, MedA.raw, MNA.raw, file="./Data/Processed Data/006.KTarray.Apac.clinical_and_array.RData")
```

####################################################################################################
##set up changeable variables
```r
expr                   = intensity.MedA
test.spots1            = apac_spots$Unique.ID[apac_spots$Unique.ID %in% apac_spots$Unique.ID[grep("std", apac_spots$Unique.ID)]==F]
test.samples           = apac_clinical$SampleID
apac1.samples          = apac_clinical$SampleID[apac_clinical$study=="Apac_X1"]
apac2.samples          = apac_clinical$SampleID[apac_clinical$study=="Apac_X2"]
apac3.samples          = apac_clinical$SampleID[apac_clinical$study=="Apac_X3"]
```

####################################################################################################
##create variable to rank antigen spots by
```r
temp <- apac_spots[apac_spots$Unique.ID %in% test.spots1,c("Unique.ID", "mean_intensity")]
temp <- temp[temp$Unique.ID %in% c(temp$Unique.ID[grep("blank", temp$Unique.ID)], temp$Unique.ID[grep("PBS", temp$Unique.ID)])==F,]
temp <- temp[order(signif(-temp$mean_intensity, digits=2)),]
temp$intensity_order <- 1:nrow(temp)

apac_spots <- merge(apac_spots, temp[c("Unique.ID", "intensity_order")], by="Unique.ID")
ag_order <- apac_spots[,c("Unique.ID","antigen","mean_intensity","intensity_order","description")]
```

####################################################################################################
##sort intensity data by sample set (apac_X1, apac_X2, apac_X3) and age:
```r
age_order <- apac_clinical[,c("sampleID","study", "Sample","age")]
age_order$random_number <- sample(1:length(age_order$Sample))
age_order <- age_order[order(age_order$study, -age_order$age, age_order$Sample),]# age_order$random_number),]
age_order$age_order2 <- 1:nrow(age_order)
age_order <- age_order[order(age_order$age_order2),]
```

####################################################################################################
##create age categories, colored by rcolorbrewer
```r
age_order$agecat <- brewer.pal(6,"RdPu")[1]
age_order$agecat[age_order$age>=3 & age_order$age<5] <- brewer.pal(6,"RdPu")[2]
age_order$agecat[age_order$age>=5 & age_order$age<7] <- brewer.pal(6,"RdPu")[3]
age_order$agecat[age_order$age>=7 & age_order$age<9] <- brewer.pal(6,"RdPu")[4]
age_order$agecat[age_order$age>=9 & age_order$age<11] <- brewer.pal(6,"RdPu")[5]
age_order$agecat[age_order$age>11] <- brewer.pal(6,"RdPu")[6]
table(age_order$agecat, age_order$age)
```

####################################################################################################
##sort subjects according to age_order
```r
expr <- intensity.MedA[rownames(intensity.MedA) %in% c(rownames(intensity.MedA)[grep("std", rownames(intensity.MedA))], 
                       rownames(intensity.MedA)[grep("PBS", rownames(intensity.MedA))], 
                       rownames(intensity.MedA)[grep("blank", rownames(intensity.MedA))])==F,]
data <- data.frame(t(expr))
data$sampleID <- rownames(data)
data <- merge(data, age_order[c("sampleID", "study", "Sample", "age", "age_order2", "agecat")], by="sampleID")
data <- data[order(data$age_order2),]
rownames(data) <- data$sampleID
data.1 <- data
```

####################################################################################################
##determine quantile breaks & color palette:
```r
temp <- data.1
min(temp[,2:(ncol(temp)-6)], na.rm=T)

quantile.range <- quantile(as.matrix(t(temp[,2:(ncol(temp)-6)])), probs = seq(0, 1, 0.01), na.rm=T)
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.7)
palette.breaks[1] <- 3.5
palette.breaks[length(palette.breaks)] <- 19
color.palette  <- brewer.pal(9,"YlGnBu")
```

####################################################################################################
##get data into the right order
```r
data <- data.1
expr_order <- data.frame(expr.index=2:(ncol(data)-5), Unique.ID=names(data)[2:(ncol(data)-5)])
ag_order2 <- merge(ag_order, expr_order, by="Unique.ID")
ag_order2 <- ag_order2[order(ag_order2$mean_intensity, decreasing=F),]
data2 <- as.matrix(data[,ag_order2$Unique.ID])
dim(data2)
#[1] 990 126
```

####################################################################################################
##create the layout
```r
tiff(filename = "./Figures/Exploratory/apac.Heatmap.tiff",
     width = 4.5, 
     height = 3.75, 
     units = "in",
     pointsize = 9,
     compression = "lzw",
     bg = "white", 
     res = 1000)
layout(matrix(c(1,1,2,3,3,4,5,5,6), 3,3, byrow = T),
       heights=c(1), widths=c(4,1))
layout.show(6)
```

####################################################################################################
##create heatmap
```r
pdf(file="./Figures/Exploratory/0006.KTarray.Apac.heatmap.pdf",width=12,height=9)
par(mfrow=c(1,2))
```
###legend for intensity
```r
par(mar=c(0.25,0.5,1,2))
plot.new()
mtext("Ab Intensity", side=1, cex=2, font=1, adj=0.5, srt=90, line=-1)
polygon(c(.4,.4,1.075),  c(.1,1,1),  border=F, col = color.palette[9])
polygon(c(.4,.4,1),     c(.1,.9,.9),  border=F, col = color.palette[8])
polygon(c(.4,.4,0.925), c(.1,.8,.8), border=F, col = color.palette[7])
polygon(c(.4,.4,0.85),  c(.1,.7,.7), border=F,  col = color.palette[6])
polygon(c(.4,.4,0.775), c(.1,.6,.6), border=F, col = color.palette[5])
polygon(c(.4,.4,0.7),   c(.1,.5,.5), border=F,  col = color.palette[4])
polygon(c(.4,.4,0.625), c(.1,.4,.4), border=F, col = color.palette[3])
polygon(c(.4,.4,0.55),  c(.1,.3,.3), border=F,   col = color.palette[2])
polygon(c(.4,.4,0.475), c(.1,.2,.2), border=F, col = color.palette[1])
```

####################################################################################################
####legend for age
```r
plot.new()
par(mar=c(0.25,0.5,1.5,2.5))
mtext("Age (years)", side=1, cex=2, font=1, adj=0.5, srt=90, line=-1)
polygon(c(.4,.4,1),    c(.1,1.05,1.05), border=F, col=brewer.pal(6,"RdPu")[6])
polygon(c(.4,.4,0.9),    c(.1,0.9,0.9), border=F, col=brewer.pal(6,"RdPu")[5])
polygon(c(.4,.4,0.8),    c(.1,0.75,0.75), border=F, col=brewer.pal(6,"RdPu")[4])
polygon(c(.4,.4,0.7), c(.1,0.6,0.6), border=F, col=brewer.pal(6,"RdPu")[3])
polygon(c(.4,.4,0.6),  c(.1,0.45,0.45), border=F, col=brewer.pal(6,"RdPu")[2])
polygon(c(.4,.4,0.5), c(.1,0.3,0.3), border=F, col=brewer.pal(6,"RdPu")[1])
text(.25, 0.95,">25", cex=1.5, font=1, xpd=T, adj=0.5)
text(.25, 0.8,"9-10", cex=1.5, font=1, xpd=T, adj=0.5)
text(.25, 0.65,"7-8", cex=1.5, font=1, xpd=T, adj=0.5)
text(.25, 0.5,"5-6", cex=1.5, font=1, xpd=T, adj=0.5)
text(.25, 0.35,"3-4", cex=1.5, font=1, xpd=T, adj=0.5)
text(.25, 0.2,"1-2", cex=1.5, font=1, xpd=T, adj=0.5)
```

####################################################################################################
##(1) create the heatmap for APAC timepoint 1
```r
par(mfrow=c(1,1))

par(mar=c(10,6,1,0.5))
image(data2[grep("X1", rownames(data2)),], col=color.palette, xlab=NA, ylab=NA, 
      xaxt='n', yaxt='n', y=1:ncol(data2), x=1:nrow(data2[grep("X1", rownames(data2)),]),
      bty='n', breaks= palette.breaks)
mtext(expression(italic(Pf)~ proteins), side=2, cex=1, font=1, adj=0.5, srt=90, line=2.5)
mtext("sorted by decreasing mean Ab intensity", side=2, cex=1, font=1, adj=0.5, srt=90, line=1)
mtext("Apac #1: Subjects, sorted by decreasing age,", side=1, cex=1, font=1, adj=0.5, line=5)
```
##########################
###draw rectangles colored according to agecats
```r
temp_data <- data[grep("X1", rownames(data)),]
for(i in 1:nrow(data2[grep("X1", rownames(data2)),])){
  rect(i-0.5, -5, i+0.5, -10,
       col=temp_data[i,"agecat"], xpd=T, border=NA)
}
```

####################################################################################################
##(2) create the heatmap for APAC timepoint 2
```r
par(mar=c(10,6,1,0.5))
image(data2[grep("X2", rownames(data2)),], col=color.palette, xlab=NA, ylab=NA, 
      xaxt='n', yaxt='n', y=1:ncol(data2), x=1:nrow(data2[grep("X2", rownames(data2)),]),
      bty='n', breaks= palette.breaks)
mtext(expression(italic(Pf)~ proteins), side=2, cex=1, font=1, adj=0.5, srt=90, line=2.5)
mtext("sorted by decreasing mean Ab intensity", side=2, cex=1, font=1, adj=0.5, srt=90, line=1)
mtext("Apac #2: Subjects, sorted by decreasing age,", side=1, cex=1, font=1, adj=0.5, line=5)
```
##########################
###draw rectangles colored according to agecats
```r
temp_data <- data[grep("X2", rownames(data)),]
for(i in 1:nrow(data2[grep("X2", rownames(data2)),])){
  rect(i-0.5, -5, i+0.5, -10,
       col=temp_data[i,"agecat"], xpd=T, border=NA)
}
```

####################################################################################################
#(3) create the heatmap for APAC timepoint 3
```r
par(mar=c(10,6,1,0.5))
image(data2[grep("X3", rownames(data2)),], col=color.palette, xlab=NA, ylab=NA, 
      xaxt='n', yaxt='n', y=1:ncol(data2), x=1:nrow(data2[grep("X3", rownames(data2)),]),
      bty='n', breaks= palette.breaks)
mtext(expression(italic(Pf)~ proteins), side=2, cex=1, font=1, adj=0.5, srt=90, line=2.5)
mtext("sorted by decreasing mean Ab intensity", side=2, cex=1, font=1, adj=0.5, srt=90, line=1)
mtext("Apac #3: Subjects, sorted by decreasing age,", side=1, cex=1, font=1, adj=0.5, line=5)
```
##########################
###draw rectangles colored according to agecats
```r
temp_data <- data[grep("X3", rownames(data)),]
for(i in 1:nrow(data2[grep("X3", rownames(data2)),])){
  rect(i-0.5, -5, i+0.5, -10,
       col=temp_data[i,"agecat"], xpd=T, border=NA)
}
dev.off()
```
