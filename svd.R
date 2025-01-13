
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
write.csv(svd.fern$v, file="fern-svec.csv",row.names=TRUE)
colnames(svd.fern$u) <- rep_len(sprintf("envsv%d",sv_ids),ncol(svd.fern$u)) #Must set row and column names prior to write.csv
rownames(svd.fern$u) <- env_vars
write.csv(svd.fern$u, file="env-svec.csv",row.names=TRUE)

df.sv <- data.frame(c(1:length(svd.fern$d)),svd.fern$d)
colnames(df.sv) <- c("index","sv")
write.csv(df.sv, file="sv.csv",row.names=FALSE)

#percents <- c("0%","20%","40%","60%","80%","100%")

mar.default <- c(5,4,4,2)+0.1 #Default margin parameter
mar.custom <- mar.default+3

pdf("cumtotalcrosscov.pdf",width=18.75,height=18.75) #pdf dimensions in inches
par(mai = c(3,3,3,2), cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
plot(cumsum(svd.fern$d^2)/sum(svd.fern$d^2),main="Cumulative percentage of cross-covariance",
     xlab="Index",ylab="",ylim=c(0,1),xlim=c(0,length((svd.fern$d))),yaxt="n",las=1)
axis(side=2,at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),las=1)
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

leftsveccoeffplot <- function(j,k,xmat,ymat){ #Plot first k coefficients of jth left singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  coeffout <- min(k,dim(svd.mats$u)[2])
  #Pick out most important components of left (environment) singular vectors
  sort.u <- sort(svd.mats$u[,j]^2,decreasing=TRUE,index.return=TRUE) #Sort according to magnitude
  sort.u$ivars <- env_vars[sort.u$ix] #Names of top envionmental factors
  topindex.u <- sort(sort.u$ix[1:coeffout]) #Sorted Indices of top coeffout most important environmental factors on jth axis
  top.u <- svd.mats$u[topindex.u,j] #Values of top coefficients
  lab.top.u <- env_vars[topindex.u] #Names following topindex.u
  par(mar=mar.custom,cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  bp.env <- barplot(height=top.u,names.arg=lab.top.u,xlim=c(-1,1),
                    main=sprintf("%s environmental axis",ordmap(j)),
                    horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Environmental factor",las=1)
  axis(side=2,at=bp.env,labels=lab.top.u,pos=-0.5,las=1) #Axis on left
}

leftsveccoeffsqplot <- function(j,k,xmat,ymat){ #Plot first k squared coefficients of jth left singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  #Pick out most important components of left (environment) singular vectors
  sort.u <- sort(svd.mats$u[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.u$ivars <- env_vars[sort.u$ix]
  coeffout <- min(k,length(sort.u$ix))
  par(mar = c(5,0,0,0)+mar.custom,cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  plot(sort.u$x[1:coeffout],ylab="Squared Coefficient",main=sprintf("%s environmental axis",ordmap(j)),
       xaxt="n",xlab="",las=3)
  axis(1, at = seq(1,coeffout,1),labels = sort.u$ivars[1:coeffout],las=3)
}

leftsveccoeffcumsumsqplot <- function(j,xmat,ymat){ #Plot cumulative contribution of squared coefficients of jth left singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  sort.u <- sort(svd.mats$u[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.u$ivars <- env_vars[sort.u$ix]
  #coeffout <- min(k,length(sort.u))
  par(mar = c(5,0,0,0)+mar.custom,cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  plot(cumsum(sort.u$x),ylim=c(0,1),ylab="Cumulative percentage sum of squares",
       main=sprintf("%s environmental axis",ordmap(j)),xaxt="n",yaxt="n",xlab="")
  #axis(1, at = c(1:length(sort.u$x)),labels = sort.u$ivars,las=2)
  axis(1, at = c(1:length(sort.u$x)),labels = sort.u$ivars,las=3)
  axis(2, at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),las=3)
  abline(h=0.9,lty=4,col="#ff0000")
  abline(h=0.8,lty=2,col="#ff0000")
}

rightsveccoeffplot <- function(j,k,xmat,ymat){ #Plot coefficients of jth right singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  coeffout <- min(k,dim(svd.mats$v[2]))
  #Pick out most important components of right singular vectors
  sort.v <- sort(svd.mats$v[,j]^2,decreasing=TRUE,index.return=TRUE) #Sort according to magnitude
  sort.v$ivars <- fern_names[sort.v$ix] #Names of top envionmental factors
  topindex.v <- sort(sort.v$ix[1:coeffout]) #Sorted Indices of top coeffout most important environmental factors on jth axis
  top.v <- svd.mats$v[topindex.v,j] #Values of top coefficients
  lab.top.v <- fern_names[topindex.v]
  par(mar=mar.custom,cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  bp.fern <- barplot(height=top.v, names.arg=lab.top.v,xlim=c(-1,1),
                    main=sprintf("%s fern species composition axis",ordmap(j)),
                    horiz=TRUE,yaxt="n",xlab="Coefficient",ylab=sprintf("Top %d species",k),las=1)
  #axis(side=2,at=bp.fern,labels=sort.v$ivars[1:coeffout],pos=-0.5,las=1) #Axis on left
  axis(side=4,at=bp.fern,labels=lab.top.v,pos=0.3,las=1) #Axis on right
}

rightsveccoeffsqplot <- function(j,k,xmat,ymat){ #Plot squared coefficients of jth right singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  coeffout <- min(k,dim(svd.mats$v[2]))
  #Pick out most important components of right singular vectors
  sort.v <- sort(svd.mats$v[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.v$ivars <- fern_names[sort.v$ix]
  #topindex.v <- sort(sort.v$ix[1:coeffout]) #Indices of top coeffout most important ferns on jth axis
  #top.v <- svd.fern$v[topindex.v,j]
  par(mar = c(5,0,0,0)+mar.custom,cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  plot(sort.v$x[1:coeffout],ylab="Squared Coefficient",main=sprintf("%s species composition axis",ordmap(j)),
       xaxt="n",xlab="",las=3)
  axis(1, at = seq(1,coeffout,1),labels = sort.v$ivars[1:coeffout],cex.lab=1.5,cex.axis=1.5,las=3)
}

rightsveccoeffcumsumsqplot <- function(j,k,xmat,ymat){ #Plot squared coefficients of jth right singular vector
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  sort.v <- sort(svd.mats$v[,j]^2,decreasing=TRUE,index.return=TRUE)
  sort.v$ivars <- fern_names[sort.v$ix]
  par(mar = c(5,0,0,0)+mar.custom,cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  plot(cumsum(sort.v$x),ylim=c(0,1),ylab="Cumulative percentage sum of squares",
       main=sprintf("%s species composition axis",ordmap(j)),yaxt="n",xlab="Index")
  #axis(1, at = c(1:length(sort.v$x)),labels = sort.v$ivars,cex.axis=1.5,cex.lab=1.5,las=1)
  axis(2, at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),las=3)
  abline(h=0.9,lty=4,col="#ff0000")
  abline(h=0.8,lty=2,col="#ff0000")
}


for(i in c(1,2)){
  pdf(sprintf("envaxis%d.pdf",i),width=18.75,height=18.75)
  leftsveccoeffplot(i,99,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("envsqcoeffs%d.pdf",i),width=18.75,height=18.75)
  leftsveccoeffsqplot(i,99,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("envcumsumsq%d.pdf",i),width=18.75,height=18.75)
  leftsveccoeffcumsumsqplot(i,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("fernaxis%d.pdf",i),width=18.75,height=18.75)
  rightsveccoeffplot(i,15,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("fernsqcoeffs%d.pdf",i),width=18.75,height=18.75)
  rightsveccoeffsqplot(i,15,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
  pdf(sprintf("ferncumsumsq%d.pdf",i),width=18.75,height=18.75)
  rightsveccoeffcumsumsqplot(i,xmat=mat.env_data,ymat=mat.bianca)
  dev.off()
}

# Note:
# Only display 15 most significant fern species for each principal component as there are 
# too many different species!
# But these only explain about 50% of the cross-covariance contributed by fern proportions.
# For 80-90%, need around 55/87 species, which is a bit messy when plotted

#Get most significant species.
#pdf("ferncumsumssq1.pdf",width=18.75,height=18.75)
#rightsveccoeffcumsumsqplot(1,xmat=mat.env_data,ymat=mat.bianca)
#dev.off()

#sort.v1 <- sort(svd.fern$v[,1]^2,decreasing=TRUE,index.return=TRUE)
#sort.v1$ivars <- fern_names[sort.v1$ix]

#Step 3: Projection plots
ftoc <- function(ftype){ #Color based on forest types
  switch(ftype, "HF"="#ff0000", "PWF"="#000080", "MDF"="#32CD32")
}
pcolors <- sapply(env_data$forest.type,ftoc)

writepvalue <- function(p){ #float input, string output
  if(p < 10^-3){ #Do we need precision?
    return("(p-value < 0.001)")
  } else {
    return(sprintf())
  }
  
}



projplot <- function(j, xmat, ymat) { #data projection plot along jth principal axes of svd, specialized for this dataset
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  yproj <- c()
  xproj <- c()
  for(i in 1:length(mat.bianca[,1])){
    yproj <- append(yproj,sum(ymat[i,]*svd.mats$v[,j])) #Right projections on vertical axis
    xproj <- append(xproj,sum(xmat[i,]*svd.mats$u[,j])) #Left projections on horizontal axis
  }
  pvalue <- cor.test(xproj,yproj)$p.value #p-value of product moment correlation coefficient
  #Convert from engineering notation to power 10 notation for presentation
  split.pmcc <- strsplit(sprintf("%.1e",pvalue),"e") #Split on the e
  val <- split.pmcc[[1]][1]
  power <- as.numeric(split.pmcc[[1]][2])
  pvalexp <- bquote(paste("(p-value = ",.(val),"\u00D7 10"^.(power),")"))#Expression object for reporting p-value
  
  plottitle <- sprintf("Projection along %s principal axes",ordmap(j))
  par(mar=mar.custom, cex=1.75, cex.axis=1.75, cex.lab=1.75, cex.main=2.5)
  plot(xproj,yproj,xlab="Environmental factors",ylab="Fern composition",main=plottitle,
       sub=pvalexp,col=pcolors,cex.sub=1.0)
  #mtext(pvalexp,side=3,cex=1.75)
  text(xproj,yproj,labels=env_data$plot,col=pcolors,cex=1.5)#LABEL PLOT IDS, THIS CAN MAKE PLOT A BIT MESSY
  legend("bottomright",legend=c("HF","PWF","MDF"),fill=c("#ff0000","#000080","#32CD32"),
         title = "Forest type",cex=1.5)
         #pch=c(1,1,1),title="Forest type",cex=1.5)
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

