
#Import fern data

env_data <- read.csv("03_Environamental_Data.csv", header=TRUE, fileEncoding="UTF-8-BOM") #Recorded environmental data
colnames(env_data) <- make.names(colnames(env_data)) #Clean up column names
#R does not like spaces in variable names
env_names <- c("site name","plot","forest type","GWC (%)","pH","N","P","% sand","OM (%)","Distance from water (m)") #Manually set names for clean presentation

#Number of individuals and number of species is already encoded in fern numbers anyways so these are somewhat redundant

bianca <- read.csv("01_Ferns_abundance_plots-no lycoph.csv",header=TRUE, fileEncoding="UTF-8-BOM") #Observed fern numbers
#Cleaning dataset
bianca <- Filter(function(x)!all(is.na(x)), bianca) #Remove empty columns. None numerics are coerced into NAs
bianca <- bianca[bianca$Family!="",] #Remove empty rows. None strings are coerced into empty strings.
bianca[is.na(bianca)] = 0 #ASSUME REMAINING EMPTY ENTRIES IN THE EXCEL FILE ARE ZEROES
colnames(bianca) <- make.names(colnames(bianca)) #Clean up column names


#get numeric columns only
num.env_data <- env_data[unlist(lapply(env_data,is.numeric))]
num.bianca <- bianca[unlist(lapply(bianca,is.numeric))]

#get names for convenience
#env_vars <- colnames(num.env_data) #This one is messy
env_vars <- env_names[unlist(lapply(env_data,is.numeric))]
fern_names <- bianca$FERN.S.scientific.name
#fern_names <- colnames(num.bianca)

# Step 1: Normalize and calculate cross covariance matrix

normalize <- function(x) {
  if(is.numeric(x)){ 
    scale(x, center=TRUE, scale=TRUE) #make sure center=TRUE to calculate deviation values, set scale=TRUE for unit variance.
  } else x
}

mat.bianca <- apply(num.bianca,1,function(x)(x/sum(x))) #Matrix of fern proportions
mat.env_data <- apply(num.env_data,2,normalize) #Matrix of normalized environmental data

cov.fern <- cov(x=mat.env_data,y=mat.bianca) #covariance matrix between normalized environmental data and fern proportions


# Step 2: Calculate SVD

svd.fern <- svd(cov.fern)

#Export SVD matrices into .csv for reference
sv_ids <- 1:ncol(svd.fern$v)
#colnames(svd.fern$v) <- paste(rep_len("fernsv",ncol(svd.fern$v)),sv_ids,sep="") #Must set row and column names prior to write.csv
colnames(svd.fern$v) <- rep_len(sprintf("fernsv%d",sv_ids),ncol(svd.fern$v))
rownames(svd.fern$v) <- fern_names
write.csv(t(svd.fern$v), file="fern-sv.csv",row.names=TRUE)
colnames(svd.fern$u) <- rep_len(sprintf("envsv%d",sv_ids),ncol(svd.fern$u)) #Must set row and column names prior to write.csv
rownames(svd.fern$u) <- env_vars
write.csv(t(svd.fern$u), file="env-sv.csv",row.names=TRUE)

df.sv <- data.frame(c(1:length(svd.fern$d)),svd.fern$d)
colnames(df.sv) <- c("index","sv")
write.csv(df.sv, file="sv.csv",row.names=FALSE)

#percents <- c("0%","20%","40%","60%","80%","100%")

mar.default <- c(5,4,4,2)+0.1 #Default margin parameter
mar.custom <- mar.default+3

pdf("cumtotalcrosscov.pdf",width=18.75,height=18.75) #pdf dimensions in inches
par(mai = c(3,3,3,2))
plot(cumsum(svd.fern$d^2)/sum(svd.fern$d^2),main="Cumulative percentage of cross-covariance",
     xlab="Index",ylab="",ylim=c(0,1),xlim=c(0,length((svd.fern$d))),yaxt="n",las=1,
     cex=3.0,cex.axis=2.5,cex.lab=2.5,cex.main=5.0)
axis(side=2,at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),las=1,cex.axis=2.5)
#axis(side=1,at=seq(0.0,7.0,1.0),labels=seq(0.0,7.0,1.0),las=1,cex.axis=3.0,pos=-0.1)
abline(h=0.9,col="#ff0000",lty=4,lwd=2.0)
abline(h=0.8,col="#ff0000",lty=2,lwd=2.0)
dev.off()

#First and second axes describe nearly 90% of variation.

#u: columns are left singular vectors (environment)
#v: columns are right singular vectors (fern composition)

ordmap <- function(j){
  #switch(j, 1="1st", 2="2nd", 3="3rd", else)
  if(j == 1){
    return("1st")
  } else if(j == 2){
    return("2nd")
  } else if(j == 3){
    return("3rd")
  } else {
    return(sprintf("%dth",j))
  }
}

leftsveccoeffplot <- function(j,xmat,ymat){ #Plot coefficients of jth left singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  par(mar=mar.custom)
  bp.env <- barplot(height=svd.fern$u[,j],names.arg=env_vars,xlim=c(-1,1),
                    main=sprintf("%s environmental axis",ordmap(j)),
                    horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Environmental factor",las=1,
                    cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0)
  axis(side=2,at=bp.env,labels=env_vars,pos=-0.5,cex.axis=2.5,las=1) #Axis on left
}

leftsveccoeffsqplot <- function(j,xmat,ymat){ #Plot squared coefficients of jth left singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  #Pick out most important components of left (environment) singular vectors
  sort.u <- sort(svd.fern$u[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.u$ivars <- env_vars[sort.u$ix]
  par(mar = c(5,0,0,0)+mar.custom)
  plot(sort.u$x,ylab="Squared Coefficient",main=sprintf("%s environmental axis",ordmap(j)),xaxt="n",xlab="",
       cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0)
  axis(1, at = c(1:length(sort.u$x)),labels = sort.u$ivars,cex.lab=1.5,cex.axis=1.5,las=2)
}

leftsveccoeffcumsumsqplot <- function(j,xmat,ymat){ #Plot squared coefficients of jth left singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  sort.u <- sort(svd.fern$u[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.u$ivars <- env_vars[sort.u$ix]
  par(mar = c(5,0,0,0)+mar.custom)
  plot(cumsum(sort.u$x),ylim=c(0,1),ylab="Cumulative percentage sum of squares",
       main=sprintf("%s environmental axis",ordmap(j)),xaxt="n",yaxt="n",xlab="",
       cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0)
  axis(1, at = c(1:length(sort.u$x)),labels = sort.u$ivars,,cex.axis=1.5,cex.lab=1.5,las=2)
  axis(2, at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),cex.axis=1.5,cex.lab=1.5,las=1)
  abline(h=0.9,lty=4,col="#ff0000")
  abline(h=0.8,lty=2,col="#ff0000")
}


rightsveccoeffplot <- function(j,xmat,ymat){ #Plot coefficients of jth right singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  bp.env <- barplot(height=svd.fern$v[,j],names.arg=env_vars,xlim=c(-1,1),
                    main=sprintf("%s environmental axis",ordmap(j)),
                    horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Environmental factor",las=1,
                    cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0)
  axis(side=2,at=bp.env,labels=env_vars,pos=-0.5,cex.axis=1.5,las=1) #Axis on left
}

rightsveccoeffsqplot <- function(j,xmat,ymat){ #Plot squared coefficients of jth right singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  #Pick out most important components of right singular vectors
  sort.v <- sort(svd.fern$v[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.v$ivars <- env_vars[sort.v$ix]
  par(mar = c(5,0,0,0)+mar.custom)
  plot(sort.v$x,ylab="Squared Coefficient",main=sprintf("%s environmental axis",ordmap(j)),xaxt="n",xlab="",
       cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0)
  axis(1, at = c(1:length(sort.v$x)),labels = sort.v$ivars,cex.lab=1.5,cex.axis=1.5,las=1)
}

rightsveccoeffcumsumsqplot <- function(j,xmat,ymat){ #Plot squared coefficients of jth right singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  sort.v <- sort(svd.fern$v[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.v$ivars <- env_vars[sort.v$ix]
  par(mar = c(5,0,0,0)+mar.custom)
  plot(cumsum(sort.v$x),ylim=c(0,1),ylab="Cumulative percentage sum of squares",
       main=sprintf("%s environmental axis",ordmap(j)),xaxt="n",yaxt="n",xlab="",
       cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0)
  axis(1, at = c(1:length(sort.v$x)),labels = sort.v$ivars,cex.axis=1.5,cex.lab=1.5,las=1)
  axis(2, at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),cex.axis=1.5,cex.lab=1.5,las=1)
  abline(h=0.9,lty=4,col="#ff0000")
  abline(h=0.8,lty=2,col="#ff0000")
}


for(i in c(1,2)){
  pdf(sprintf("envaxis%d.pdf",i),width=18.75,height=18.75)
  leftsveccoeffplot(i,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("envsqcoeffs%d.pdf",i),width=18.75,height=18.75)
  leftsveccoeffsqplot(i,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("envcumsumsq%d.pdf",i),width=18.75,height=18.75)
  leftsveccoeffcumsumsqplot(i,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
}

#Note: Fern compositions omitted here because a relatively large proportion of species is needed to explain each principal axis.

#Step 3: Projection plots
ftoc <- function(ftype){ #Color based on forest types
  switch(ftype, "HF"="#ff0000", "PWF"="#000080", "MDF"="#32CD32")
}
pcolors <- sapply(env_data$forest.type,ftoc)

projplot <- function(j, xmat, ymat) { #data projection plot along jth principal axes of svd, specialized for this dataset
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  yproj <- c()
  xproj <- c()
  for(i in 1:length(mat.bianca[,1])){
    yproj <- append(yproj,sum(ymat[i,]*svd.mats$v[,j])) #Right projections on vertical axis
    xproj <- append(xproj,sum(xmat[i,]*svd.mats$u[,j])) #Left projections on horizontal axis
  }
  plottitle <- sprintf("Projection along %s principal axes",ordmap(j))
  par(mar=mar.custom)
  plot(xproj,yproj,xlab="Environmental factors",ylab="Fern composition",main=plottitle,col=pcolors,
       cex=2.5,cex.axis=2.5,cex.lab=3.0,cex.main=5.0) 
  text(xproj,yproj,labels=env_data$plot,col=pcolors,cex=2.5)#LABEL PLOT IDS, THIS CAN MAKE PLOT A BIT MESSY
  legend("bottomright",legend=c("HF","PWF","MDF"),col=c("#ff0000","#000080","#32CD32"),
         pch=c(1,1,1),title="Forest type",cex=2.5)
  abline(h=0)
  abline(v=0)
} #BUG: Plot IDs are not labelled FIXED

pdf("axproj1.pdf",width=18.75,height=18.75)
projplot(1, xmat=mat.env_data, ymat=mat.bianca)
dev.off()
pdf("axproj2.pdf",width=18.75,height=18.75)
projplot(2, xmat=mat.env_data, ymat=mat.bianca)
dev.off()
#projplot(3, xmat=mat.env_data, ymat=mat.bianca) #Projection plot looks too mixed to visualize any meaningful difference.

