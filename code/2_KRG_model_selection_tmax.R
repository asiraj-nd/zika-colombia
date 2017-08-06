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
temp = read.csv("../data/station_tmax.csv", header=T)
numdays<-dim(temp)[2]-5 # droping label columns 
locator<-temp[,2:5]  # label columns
loct  = read.csv("../data/station_serial.csv")[,2]

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
    fdata = paste("../data/wclim/tmax",mon,"_col.csv",sep="")
    tmp <- as.vector(unlist(read.csv(fdata, header=T)))
  }  
  tmpthis = tmp[loct]
  allddwctemp<- cbind(allddwctemp,tmpthis)
  prvmonYr = mon
}

###########
n.cores <- 30
n.lines <- 30

seqP <- seq(1,n.lines,n.lines/n.cores)
xx <- as.numeric(Sys.getenv("SGE_TASK_ID"))
row.num <- seqP[xx]

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
                          locator[-drop.station,4],allddwctemp[-drop.station,coler-5],allddtmp[-drop.station,coler-5])) 
    myd = myd [which(!is.na(myd[,3])),]
    names(myd)<- c("lon","lat","temp", "dist", "vqx","wqx","nqx")
    head(myd)
    row.names(myd) = NULL
    if (is.na(temp[drop.station,coler])) {
      finalout = NA
    } else {
      
      mygeod=  as.geodata(myd[,c(1,2,3,5,6,7)], coords.col=1:2, data.col=3, covar.col=c(1:2,4:6))
      maxd = 1.5*diff(range(myd$lat))
      var2 <- variog(mygeod, option="bin",
                     trend=~lon+lat+vqx+wqx+nqx,
                     bin.cloud="TRUE", max.dist=maxd)
      
      mxr = mean(var2$v[order(-var2$v)][1:3]) #mean of the top three
      fit3 = try(variofit(var2, cov.model="gauss", ini.cov.pars=c(mxr,8), 
                          fix.nugget= FALSE, nugget=10, wei="cressie"), silent=T)
      
      if(is(fit3,"try-error"))  {
        mxr =0.0005
        nng = 0.0
      } else { nng = fit3$nugget
      if (nng<0) nng=0 }
      monYr = allmn[allddmon[coler-5,4]]
      mon = allddmon[coler-5,2]
      
      if (prvmonYr!=monYr) {  
        fdata = paste("../data/wclim/tmax",mon,"_col.bil",sep="")
        tmpwclim <- as.matrix(raster(fdata),ncol=nx,nrow=ny)
        min(tmpwclim,na.rm=T)
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
      #image.plot( as.image(kc$predict,mygrid.data[,c("lon","lat")],
      #                    nrow=480, ncol=432))
      #quilt.plot(x,y)            
      #image(kc, loc = tolocs, col=gray(seq(1,0.1,l=30)), xlab="Coord X", ylab="Coord Y")
      #max(kc$predict)
      
      kkc= matrix(kc$predict,nrow=480)
      finalout = reverse_grid_match(kkc)
      
      # finalout[which(finalout<0)] = 0 #for rainfall
      #zq<- structure(.Data=as.vector(unlist(t(finalout[480:1,]))), .Dim=c(nx,ny))
      #zq<-list(x=xq, y=yq, z=zq)
      #image(zq)
    }      
    thisstation <- c(thisstation,finalout[loct[drop.station]])  # extract prediction at dropped station
    prvmonYr = mon
  }
  txt<-paste("../generated/interpol_tmax_wo_stn",drop.station,"_allsurf.csv",sep="")
  write.csv(thisstation,txt,row.names=F, quote=F)
}


  ##### generate statistic
  obspred<- NULL
  for (drop.station in 1:30){
    dd = as.vector(unlist(read.csv(paste("../data/climinterpolated/interpol_tmax_wo_stn",drop.station,"_allsurf.csv",sep=""))))
    obspred = rbind(obspred,cbind(as.vector(unlist((temp[drop.station,-(1:5)]))),round(dd,1)))
  }  

  obspred[,1] = (obspred[,1] -32)* 5/9  #convert temp to oC
  eee<- obspred[apply(obspred,1,function(tr){!any(is.na(tr))}),]  # drop NA

  mae(eee[,1], eee[,2])
  cor(eee[,1], eee[,2])
  cv(eee[,1], eee[,2])

