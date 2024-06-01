BiocManager::install("GEOquery", force = TRUE)

library(GEOquery)

if(!file.exists("geo_downloads")) 
  dir.create("geo_downloads")

my.gse <- "GSE12195"  

if(!file.exists(paste0("./geo_downloads/",my.gse)))
  getGEOSuppFiles(my.gse, makeDirectory=T, baseDir="geo_downloads")


my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="./geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)

my.geo.gse

my.geo.gse <- my.geo.gse[[1]]


my.geo.gse

untar(paste0("geo_downloads/",my.gse,"/",my.gse,"_RAW.tar"), exdir=paste0("geo_downloads/",my.gse,"/CEL/"))
my.cels <- list.files(paste0("geo_downloads/",my.gse,"/CEL/"))

my.cels

#Preparing the Phenodata
my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
head(my.pdata)
dim(my.pdata)
colnames(my.pdata)
table(rownames(my.pdata) == my.cels) # it returned false, it is normally

head(my.cels)
head(rownames(my.pdata))

temp.rownames <- rownames(my.pdata)
my.cels <- gsub(".CEL.gz", "", my.cels)
temp.rownames
my.cels
filtered.temp.rownames <- temp.rownames[temp.rownames   %in% my.cels]

filtered.temp.rownames

my.cels
table(filtered.temp.rownames == my.cels)

# Filter dataframe 'my.pdata' for rows matching 'filtered.temp.rownames'
# Update to retrieve relevant rows from original dataframe instead of NA for non-matching rows
matched_rows <- match(filtered.temp.rownames, rownames(my.pdata))
my.pdata <- my.pdata[matched_rows, ]

temp.rownames2 <- paste(rownames(my.pdata), ".CEL.gz", sep="")
rownames(my.pdata) <- temp.rownames2
# removing NA
my.pdata <- my.pdata[!is.na(matched_rows), ]

# writing new row names and data to control
rownames(my.pdata)
head(my.pdata)


# writing new row names and data to control
rownames(my.pdata)

rm(temp.rownames)
rm(filtered.temp.rownames)
table(rownames(my.pdata) == my.cels)

write.table(my.pdata, file=paste0("geo_downloads/",my.gse,"/CEL/",my.gse,"_PhenoData.txt"), sep="\t", quote=F)

#Reading the CEL Files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


# checking affy package that it is installed or not
if (!requireNamespace("affy", quietly = TRUE)) {
  BiocManager::install("affy")
} else {
  message("affy was installed.")
}

library(affy)
cel.path <- paste0("./geo_downloads/",my.gse,"/CEL")
my.affy <- ReadAffy(celfile.path=cel.path, phenoData=paste(cel.path, paste0(my.gse,"_PhenoData.txt"), sep="/"))
show(my.affy)
exprs(my.affy)
head(exprs(my.affy))

colnames(pData(my.affy))
pData(my.affy)
pData(my.affy)$title
pData(my.affy)$description

#Function that converts the Affy object into an rma object and normalizes it
my.rma <- rma(my.affy, normalize=T, background=T)  #quantile normalizasyonu.
head(exprs(my.rma))

#Annotation#
my.rma@annotation
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db) #BiocManager::install("hgu133plus2.db")
library(annotate)
library(R2HTML)

ID <- featureNames(my.rma)
Symbol<-getSYMBOL(ID,"hgu133plus2.db")
sym <- as.data.frame(Symbol)

data <- as.data.frame(exprs(my.rma)) 
data <- cbind(sym,data)  #sC<tun ekleme

i <- which(is.na(data$Symbol) == TRUE)
data<-data[-c(i),]

rownames(data) <- data[,1] 

#library(data.table) allows you to perform more operations than data.table, data.frame.
X <- data.table::as.data.table(data)
final_data <- X[,lapply(.SD,mean),"Symbol"]
final_data <- as.data.frame(final_data)
rownames(final_data) <- final_data[,1] 
final_data <- final_data[,-c(1)]

saveRDS(final_data,"geo_downloads/GSE12195/GSE12195_raw.RDS")

final_data  = t(final_data) 
metadata = pData(my.affy)

#We checked that the samples in both datasets are in the same order.
table (rownames(final_data) == rownames(metadata))

final_data = as.data.frame(final_data)
final_data$stage = metadata$title
#final_data$stage = as.numeric(final_data$stage) ##If we want the stage values ​​to be numbers, not chars, we change them like this.

write.csv(final_data,file="geo_downloads/GSE12195/GSE12195.csv")

# taking row names
row_names <- rownames(metadata)

# taking final_data data
stage_values <- final_data$stage

# saving data
write.csv(stage_values, "geo_downloads/GSE12195/GSE12195_label.csv", row.names = FALSE)

