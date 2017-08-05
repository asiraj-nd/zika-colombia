####################################################
generate.Header<-function(nx,ny,minlon,maxlat,txt,cellWidth) {
  A1=c("ENVI description", "=", "{File Resize Result, x resize factor: 1.000000, y resize factor: 1.000000. [Sun Nov 05 11:08:48 2012]}")
  A2=c("samples", "=",nx)
  A3=c("lines", "=",ny)
  A4=c("bands",  "=",1)
  A5=c("header offset",  "=",0)
  A6=c("file type", "=","ENVI Standard")
  A7=c("data type", "=",4)
  A8=c("interleave", "=","bsq")
  A9=c("sensor type", "=","Unknown")
  A10=c("byte order ", "=",0)
  A11=c("x start", "=",1)
  A12=c("y start", "=",1)
  A13=c("map info", "=",paste("{Geographic Lat/Lon, 1.0000, 1.0000,", minlon, ",",maxlat,",", cellWidth,",", cellWidth,", WGS-84, units=Meters}",sep=""))
  A14=c("wavelength units", "=","Unknown")
  A15=c("data ignore value", "=",-9.99900000e+003)
  A16=c("band names", "=","{Resize (Band 1:tmin_1.bil)}")
  write.table(rbind(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16), txt, sep=" ",col.names=F,row.names=F, quote=F)
  return(0)
}

####################################################
grid_match<- function(dem, nx=432) {
  vq<-dem[480:1,]
  zq<- structure(.Data=as.vector(unlist(t(vq))), .Dim=c(nx,480))
  return(zq)
}

### 2.5' raster output (colombia)

colo2_5<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-432
  ny<-480
  xllCorner<- -83
  maxlat<-15.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.04166667
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

reverse_grid_match<- function(fhat, ny=432) {
  finalout<-structure(fhat,.Dim=c(ny,480))
  finalout<- t(finalout)[480:1,]
  return(finalout)
}

#####
resample<- function(mxd,scl) {
  rw = dim(mxd)[1]
  cl = dim(mxd)[2]
  nwrw= rw*scl
  nwcl= cl*scl
  ecogridmx<- matrix(NA,nrow=rw, ncol=nwcl)
  for (jj in 0:(cl-1))
    ecogridmx[,(1:scl)+jj*scl]<- rep(mxd[,jj+1],scl)
  ecogridxmx<- matrix(NA,nwrw,nwcl)
  for (ii in 0:(rw-1))
    ecogridxmx[(1:scl)+ii*scl,]<- matrix(rep(ecogridmx[ii+1,],scl),scl,nwcl,byrow=T)
  return(ecogridxmx)
}  

rmse<-function(a,b)
  sqrt(sum(((a-b)^2)/length(a)))

cv<-function(a,b)
  rmse(a,b)/mean(a)

mae<-function(a,b)
  sum(abs(a-b))/ length(a)

##############################################################
generate.Header.sinu<-function(nx,ny,minlon,maxlat,txt,cellWidth) {
  A1=c("ENVI description", "=", "{File Resize Result, x resize factor: 1.000000, y resize factor: 1.000000. [Sun Nov 05 11:08:48 2012]}")
  A2=c("samples", "=",nx)
  A3=c("lines", "=",ny)
  A4=c("bands",  "=",1)
  A5=c("header offset",  "=",0)
  A6=c("file type", "=","ENVI Standard")
  A7=c("data type", "=",4)
  A8=c("interleave", "=","bsq")
  A9=c("sensor type", "=","Unknown")
  A10=c("byte order ", "=",0)
  A11=c("x start", "=",1)
  A12=c("y start", "=",1)
  A13=c("map info", "=",paste("{Sinusoidal, 1, 1,", minlon, ",",maxlat,",", cellWidth,",", cellWidth,",WSG-84 , units=Meters}",sep=""))
  A14=c("wavelength units", "=","Unknown")
  A15=c("data ignore value", "=",-9.99900000e+003)
  A16=c("band names", "=","{Resize (Band 1:tmin_1.bil)}")
  A17=c("projection info", "=",'{PROJCS["Sinusoidal",GEOGCS["GCS_Undefined",DATUM["Undefined",SPHEROID["User_Defined_Spheroid",6371007.181,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Sinusoidal"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],UNIT["Meter",1.0]]}')
  write.table(rbind(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17), txt, sep=" ",col.names=F,row.names=F, quote=F)
  return(0)
}

###############
modis.ndvi.mozaic<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-2400
  ny<-3600
  xllCorner<- -8895604.1581319998949766
  maxlat<-2223901.0393329998478293
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-926.6254331383334
  generate.Header.sinu(nx,ny,minlon,maxlat,txt,cellWidth)
}

###### aggregate to lower resolution
glbcoarser<- function(gmx,rt,aggrg) {
  gx<-dim(gmx)/rt
  shrinkmx<-matrix(NA,gx[1],gx[2])
  for (i in 0:(gx[1]-1))
    for (j in 0:(gx[2]-1)) { 
      if (aggrg==1 & any(!is.na(gmx[(1:rt)+i*rt,(1:rt)+j*rt]))) shrinkmx[i+1, j+1]<- mean(gmx[(1:rt)+i*rt,(1:rt)+j*rt],na.rm=T)
      if (aggrg==2 & any(!is.na(gmx[(1:rt)+i*rt,(1:rt)+j*rt]))) shrinkmx[i+1, j+1]<- sum (gmx[(1:rt)+i*rt,(1:rt)+j*rt],na.rm=T)
    }
  return(shrinkmx)
}

######
modis.ndvi.resample<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-438
  ny<-480
  xllCorner<-  -9229157.0434563197195530
  maxlat<-1667925.7796497857198119
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-926.6254331383334  
  generate.Header.sinu.reproj(nx,ny,minlon,maxlat,txt,cellWidth)
}

##### 1 degree raster output (global)
output1deg<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-360
  ny<-180
  maxlat<-90.0
  minlon<- -180
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-1
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

### 2.5' raster output (global)
output2_5<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-8640
  ny<-4320
  xllCorner<--180
  maxlat<-90.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.041666667
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}


### resample 1 deg to 2.5' resolution
resample_2_5<- function(mxd) {
  ecogridmx<- matrix(NA,nrow=180, ncol=8640)
  for (jj in 0:359)
    ecogridmx[,(1:24)+jj*24]<- rep(mxd[,jj+1],24)
  ecogridxmx<- matrix(NA,4320,8640)
  for (ii in 0:179)
    ecogridxmx[(1:24)+ii*24,]<- matrix(rep(ecogridmx[ii+1,],24),24,8640,byrow=T)
  return(ecogridxmx)
}  

#### raster output 0.5'
colo0_5<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-2160
  ny<-2400
  xllCorner<- -83
  maxlat<-15.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.008333
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

#### raster output 0.05'
colo0_05<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-21600
  ny<-24000
  xllCorner<- -83
  maxlat<-15.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.0008333
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

#### raster output 0.25'
colo0_25<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-4320
  ny<-4800
  xllCorner<- -83
  maxlat<-15.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidth<-0.004166667
  generate.Header(nx,ny,minlon,maxlat,txt,cellWidth)
}

######## linear interpolation
interpolateCols <- function(vec, nSubsteps) {
  matInterp <- rep(0,length(vec)+(length(vec)-1)*nSubsteps)
  for(i in 1:(length(vec)-1)) {
    interpVals <- approx(vec[i:(i+1)], n=(nSubsteps+2))
    lilloc<-i+nSubsteps*(i-1)
    matInterp[lilloc:(lilloc+nSubsteps+1)] <- interpVals$y
  }
  return(matInterp)
}

####### output nesdis rainfall raster
SA_nesdis<- function (mx,outdf1) {
  outdf = paste(outdf1,".bil",sep="")
  writeBin(as.numeric(as.vector(t(mx))), outdf, size = 4)
  print(outdf)
  nx<-1359
  ny<-1613
  xllCorner<--82
  maxlat<-13.0
  minlon<-xllCorner
  txt<-paste(outdf1,".hdr",sep="")
  cellWidthx<-0.0382513
  cellWidthy<-0.0359477
  generate.Header.uneven(nx,ny,minlon,maxlat,txt,cellWidthx,cellWidthy)
}
