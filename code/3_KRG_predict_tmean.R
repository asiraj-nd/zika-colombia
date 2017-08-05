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

### read station temperature 
temp = read.csv("../data/station_tmean.csv", header=T)
numdays<-dim(temp)[2]-5
locator<-temp[,2:5]

##### vector matching data 
Monames= c("Jan","Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
allmn= paste(Monames,"_",rep(2011:2016,12)[order(rep(2011:2016,12))], sep="")[1:69]
allddmon = read.csv("../data/date_mon.csv") 

##### reading monthly worldclim data
allddwctemp = NULL
prvmonYr = 0

for (i in 1:(dim(allddmon)[1])){
  mon = allddmon[i,2]
  if (prvmonYr!=mon) {  
    fdata = paste("../data/wclim/tmean",mon,"_col.csv",sep="")
    tmp <- as.vector(unlist(read.csv(fdata, header=T)))
  }  
  tmpthis = tmp[loct]
  allddwctemp<- cbind(allddwctemp,tmpthis)
  prvmonYr = mon
}

for (coler in 5+(1:numdays)){
    myd<-data.frame(cbind(locator[,2:1],as.vector(unlist(temp[,coler]-32)*5/9),as.vector(unlist(temp[,4])), 
                          locator[,4],allddwctemp[,coler-5])) 
    myd = myd [which(!is.na(myd[,3])),]
    names(myd)<- c("lon","lat","temp", "dist", "vqx","wqx")
    row.names(myd) = NULL

    mygeod=  as.geodata(myd[,c(1,2,3,5,6)], coords.col=1:2, data.col=3, covar.col=c(1:2,4:5))
    maxd = 1.5*diff(range(myd$lat))
    var2 <- variog(mygeod, option="bin",
                   trend=~lon+lat+vqx+wqx,
                   bin.cloud="TRUE", max.dist=maxd)

    mxr = mean(var2$v[order(-var2$v)][1:3]) #mean of the top three
    fit3 = try(variofit(var2, cov.model="gauss", ini.cov.pars=c(mxr,8), 
                        fix.nugget= FALSE, nugget=10, wei="cressie"), silent=T)
    
    if(is(fit3,"try-error"))  {
      mxr =0.0005
      nng = 0.0
    } else { nng = fit3$nugget
    if (nng<0) nng=0 }
    
    mon = allddmon[coler-5,2]
    if (prvmonYr!=mon) {  
      fdata = paste("../data/wclim/tmean",mon,"_col.csv",sep="")
      tmp <- as.vector(unlist(read.csv(fdata, header=T)))
      tmpwclim <- matrix(as.vector(unlist(tmp)),nrow=480)
    }  
    tmpwclim[is.na(tmpwclim)]<- 0
    wq<- grid_match(tmpwclim)  
    
    mygrid.data = cbind(1,tolocs,vqx=c(vq),wqx=c(wq))
    grid.geod = as.geodata(mygrid.data, coords.col=2:3, data.col=1, covar.col=2:5)
    
    kc = krige.conv(geodata=mygeod, locations=grid.geod$coords, 
                            krige=krige.control(type.krige="OK", 
                              trend.d=trend.spatial(var2$trend,mygeod),
                              trend.l=trend.spatial(~lon+lat+vqx+wqx, grid.geod),
                              cov.pars=c(mxr,8), nug=nng))

      kkc= matrix(kc$predict,nrow=480)
      finalout = reverse_grid_match(kkc)
      #zq<- structure(.Data=as.vector(unlist(t(finalout[480:1,]))), .Dim=c(nx,ny))
      #zq<-list(x=xq, y=yq, z=zq)
      #image(zq)
    range(finalout, na.rm=T)
    colo2_5(finalout,txt)
    prvmonYr = mon
}

  ### NA to ocean
  txt<-paste("../generated/tmean_day_",104,sep="")
  tmeanx <- as.matrix(raster(paste(txt,".bil",sep="")),ncol=nx,nrow=ny)
  ocean = which(tmeanx < -16, arr.ind=T)
  
  for (i in 1:1096) {
    txt<-paste("../generated/tmean_day_",i,sep="")
    tmean <- as.matrix(raster(paste(txt,".bil",sep="")),ncol=nx,nrow=ny)
    tmean[ocean]<- NA
    colo2_5(tmean,txt)
  }

