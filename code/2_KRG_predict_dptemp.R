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
temp = read.csv("../data/station_dptemp.csv", header=T)
numdays<-dim(temp)[2]-5 # droping label columns 
locator<-temp[,2:5]  # label columns

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
      txt<-paste("../generated/dptemp_day_",coler-5,sep="")
      colo2_5(finalout,txt)
}


  ### NA to ocean
  txt<-paste("../generated/tmean_day_",104,sep="")
  tmeanx <- as.matrix(raster(paste(txt,".bil",sep="")),ncol=nx,nrow=ny)
  ocean = which(tmeanx < -16, arr.ind=T)

  for (i in 1:1096) {
    txt<-paste("../generated/dptemp_day_",i,sep="")
    tmean <- as.matrix(raster(paste(txt,".bil",sep="")),ncol=nx,nrow=ny)
    tmean[ocean]<- NA
    colo2_5(tmean,txt)
  }


####### correct erroneous tmin values (where tmax<tmin)

#j=353
for (j in 967:1096){
  txt<-paste("../generated/tmean_day_",j,sep="")
  fdata = paste(txt,".bil",sep="")
  dd1 <- as.matrix(raster(fdata),ncol=nx,nrow=ny)
  txt<-paste("../generated/dptemp_day_",j,sep="")
  fdata = paste(txt,".bil",sep="")
  dd2 <- as.matrix(raster(fdata),ncol=nx,nrow=ny)
  negdif = which(dd2>dd1)
  if (length(negdif)>0) {
    dd2[negdif] = dd1[negdif]
    colo2_5(dd2,txt)
    print (c("yes",j))
  }
}


#### generate relative humidity
#j=353
for (j in 967:1096){
  txt<-paste("../generated/tmean_day_",j,sep="")
  fdata = paste(txt,".bil",sep="")
  dd1 <- as.matrix(raster(fdata),ncol=nx,nrow=ny)
  txt<-paste("../generated/dptemp_day_",j,sep="")
  fdata = paste(txt,".bil",sep="")
  dd2 <- as.matrix(raster(fdata),ncol=nx,nrow=ny)
  relh= exp(-17.271*237.7*(dd1-dd2)/(237.7+dd2)/(237.7+dd1))
  txt<-paste("../generated/relhum_day_",j,sep="")
  colo2_5(relh,txt)
  print (c("yes",j))
}


#### validate relative humidity
#j=353
allrh= NULL
for (j in 1:966){
  txt<-paste("../generated/relhum_day_",j,sep="")
  fdata = paste(txt,".bil",sep="")
  dd1 <- as.matrix(raster(fdata),ncol=nx,nrow=ny)
  allrh = cbind(allrh,dd1[loct])
  print (c("yes",j))
}

temp2 = read.csv("../data/station_humidity.csv", header=T)
temp2<- temp2[-c(21,26),] 

length(temp2)/length(dd1)

length(as.vector(unlist(temp[,-(1:5)])))
length(as.vector(unlist(allrh*100)))


ddd = cbind(as.vector(unlist(allrh*100)),as.vector(unlist(temp[,-(1:5)])))
head(ddd)
eee<- ddd[apply(ddd,1,function(tr){!any(is.na(tr))}),] 

mae(eee[,1], eee[,2])
cor(eee[,1], eee[,2])
cv(eee[,1], eee[,2])

plot(eee[,1], eee[,2],pch=19 ,  cex=.5)
abline(lm(eee[,2]~eee[,1]), col='red')
abline(a=0,b=1, col='black',lty=2)

#incl=loct[which(!is.na(as.vector(unlist(temp[,coler]))))]
#dded = incl[which(incl!=loct[drop.station])]
#layout(1)
#plot(cbind(finalout[dded],exp(tmp[,3]))-1)
#abline(a=0,b=1)
#dem2.5[dded]
