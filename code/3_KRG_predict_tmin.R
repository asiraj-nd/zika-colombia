set.seed(123)
library(geoR)
library(fields)

setwd("/afs/crc.nd.edu/user/a/asiraj/rapid/code/")
source("0_numfunctions.R")
dem2.5 = read.csv("../data/dem2_5.csv", header=T)
vq=grid_match(dem2.5)

#### set up grid space
xq<-seq(-83,-65,0.0416667)[1:432]
yq<-seq(15,-5,-0.0416667)[1:480]
yq<-yq[order(yq)]
grid.list<- list( x=xq, y= yq)
bigX<- make.surface.grid( grid.list)
cellWidth = 0.0416667
tolocs  =  expand.grid(lon=xq,lat=yq)
tolocs = as.data.frame(tolocs)
nx=length(xq)
ny=length(yq)
mygrid.data = cbind(1,tolocs)
grid.geod = as.geodata(mygrid.data, coords.col=2:3, data.col=1, covar.col=2:3)

### read station data
temp = read.csv("../data/station_tmin.csv", header=T)
numdays<-dim(temp)[2]-5 # droping label columns 
locator<-temp[,2:5]  # label columns
loct  = read.csv("../data/station_serial.csv")[,2]

##### vector matching data 
Monames= c("Jan","Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
allmn= paste(Monames,"_",rep(2011:2016,12)[order(rep(2011:2016,12))], sep="")[1:69]
allddmon = read.csv("../data/date_mon.csv") 

##### reading daily tmean data
allddtmp = NULL
for (i in 1:(dim(allddmon)[1])){
  fdata = paste("../gen2/tmean_day_",i,".csv",sep="")
  tmp <- as.vector(unlist(read.csv(fdata, header=T)))
  tmpthis = tmp[loct]
  allddtmp<- cbind(allddtmp,tmpthis)
}
allddtmp

###########
n.cores <- 30
n.lines <- 30

seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx]

allstation <- NULL
for (coler in 5+(967:numdays)){
  myd<-data.frame(cbind(locator[,2:1],as.vector(unlist(temp[,coler]-32)*5/9),as.vector(unlist(temp[,4])), 
                        locator[,4],allddtmp[,coler-5])) 
  myd = myd [which(!is.na(myd[,3])),]
  names(myd)<- c("lon","lat","temp", "dist", "vqx","nqx")
  head(myd)
  row.names(myd) = NULL
  
  mygeod=  as.geodata(myd[,c(1,2,3,5,6)], coords.col=1:2, data.col=3, covar.col=c(1:2,4:5))
  maxd = 1.5*diff(range(myd$lat))
  ## stockastic model (non- repeatable)
  var2 <- variog(mygeod, option="bin",
                 trend=~lon+lat+vqx+nqx,
                 bin.cloud="TRUE", max.dist=maxd)
  
  mxr = mean(var2$v[order(-var2$v)][1:3]) #mean of the top three
  fit3 = try(variofit(var2, cov.model="gauss", ini.cov.pars=c(mxr,8), 
                      fix.nugget= FALSE, nugget=10, wei="cressie"), silent=T)
  
  if(is(fit3,"try-error"))  {
    mxr =0.0005
    nng = 0.0
  } else { nng = fit3$nugget
  if (nng<0) nng=0 }
  
  fdata = paste("../generated/tmean_day_",coler-5,".bil",sep="")
  tmp <- as.vector(unlist(as.matrix(raster(fdata),ncol=nx,nrow=ny)))
  tmpnoaa <- matrix(as.vector(unlist(tmp)),nrow=480)
  tmpnoaa[is.na(tmpnoaa)]<- 0
  nq<- grid_match(tmpnoaa)  
  
  mygrid.data = cbind(1,tolocs,vqx=c(vq),nqx=c(nq))
  grid.geod = as.geodata(mygrid.data, coords.col=2:3, data.col=1, covar.col=2:5)
  
  kc = krige.conv(geodata=mygeod, locations=grid.geod$coords, 
                  krige=krige.control(type.krige="OK", 
                                      trend.d=trend.spatial(var2$trend,mygeod),
                                      trend.l=trend.spatial(~lon+lat+vqx+nqx,grid.geod),
                                      cov.pars=c(mxr,8), nug=nng))
  kkc= matrix(kc$predict,nrow=480)
  finalout = reverse_grid_match(kkc)
  # finalout[which(finalout<0)] = 0 #for rainfall
  #zq<- structure(.Data=as.vector(unlist(t(finalout[480:1,]))), .Dim=c(nx,ny))
  #zq<-list(x=xq, y=yq, z=zq)
  #image(zq)
  txt<-paste("../generated/tmin_day_",coler-5,sep="")
  colo2_5(finalout,txt)
}


### NA to ocean
txt<-paste("../generated/tmean_day_",104,sep="")
tmeanx <- as.matrix(raster(paste(txt,".bil",sep="")),ncol=nx,nrow=ny)
ocean = which(tmeanx < -16, arr.ind=T)

for (i in 1:1096) {
  txt<-paste("../generated/tmin_day_",i,sep="")
  tmean <- as.matrix(raster(paste(txt,".bil",sep="")),ncol=nx,nrow=ny)
  tmean[ocean]<- NA
  colo2_5(tmean,txt)
}


