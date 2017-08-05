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


allrng=NULL
###########
n.cores <- 30
n.lines <- 30

seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx]

allstation <- NULL
for (drop.station in row.num:row.num) {
  thisstation<- NULL
  prvmonYr = 0
  pars=NULL
  
  for (coler in 5+(1:numdays)){
    myd<-data.frame(cbind(locator[-drop.station,2:1],as.vector(unlist(temp[-drop.station,coler]-32)*5/9),as.vector(unlist(temp[-drop.station,4])), 
                          locator[-drop.station,4],allddtmp[-drop.station,coler-5])) 
    myd = myd [which(!is.na(myd[,3])),]
    names(myd)<- c("lon","lat","temp", "dist", "vqx","nqx")
    myd[,6] = myd[,6]
    head(myd)
    row.names(myd) = NULL
    
    if (is.na(temp[drop.station,coler])) {
      finalout = NA
    } else {
      
      mygeod=  as.geodata(myd[,c(1,2,3,5,6)], coords.col=1:2, data.col=3, covar.col=1:5)
      maxd = 1.5*diff(range(myd$lat))
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
      
      fdata = paste("../gen2/tmean_day_",coler-5,".csv",sep="")
      tmp <- as.vector(unlist(read.csv(fdata, header=T)))
      tmpnoaa <- matrix(as.vector(unlist(tmp)),nrow=480)
      tmpnoaa[is.na(tmpnoaa)]<- 0
      nq<- grid_match(tmpnoaa)  

      mygrid.data = cbind(1,tolocs,vqx=c(vq),nqx=c(nq))
      grid.geod = as.geodata(mygrid.data, coords.col=2:3, data.col=1, covar.col=2:5)
      
      kc = krige.conv(geodata=mygeod, locations=grid.geod$coords, 
                      krige=krige.control(type.krige="OK", 
                                          trend.d=trend.spatial(var2$trend,mygeod),
                                          trend.l=trend.spatial(~lon+lat+vqx+nqx, grid.geod),
                                          cov.pars=c(mxr,8), nug=nng))
      kkc= matrix(kc$predict,nrow=480)
      finalout = reverse_grid_match(kkc)
      finalout[which(finalout<0)] = 0
      #zq<- structure(.Data=as.vector(unlist(t(finalout[480:1,]))), .Dim=c(nx,ny))
      #zq<-list(x=xq, y=yq, z=zq)
      #image(zq)
    }      
    thisstation <- c(thisstation,finalout[loct[drop.station]])  # extract prediction at dropped station
  }
  txt<-paste("../generated/interpol_dptemp_wo_stn",drop.station,"_allsurf.csv",sep="")
  write.csv(thisstation,txt,row.names=F, quote=F)
}


##### generate statistic
obspred<- NULL
for (drop.station in 1:30){
  dd = as.vector(unlist(read.csv(paste("../data/climinterpolated/interpol_dptemp_wo_stn",drop.station,"_allsurf.csv",sep=""))))
  obspred = rbind(obspred,cbind(as.vector(unlist((temp[drop.station,-(1:5)]))),round(dd,1)))
}  

obspred[,1] = (obspred[,1] -32)* 5/9  #convert temp to oC
eee<- obspred[apply(obspred,1,function(tr){!any(is.na(tr))}),]  # drop NA

mae(eee[,1], eee[,2])
cor(eee[,1], eee[,2])
cv(eee[,1], eee[,2])

