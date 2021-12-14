library(lattice);
library(mapproj);
library(RNetCDF);
library(ncdf.tools);
#####################################
open_dat <- open.nc("/Users/ccfer/Downloads/air.2m.gauss.2021.nc");
global_dat <- read.nc(open_dat); 
dat <- global_dat$air;   # this is a (nlon x nlat x ntime) array with the observations
#####################################
#####################################
lat = global_dat$lat;
lon = global_dat$lon;
for(j in 1:length(lon)){ if(lon[j]>180){lon[j] = lon[j]-360;} }

coord = expand.grid(lon,lat);
#####################################

dat1 = array(dim=c(length(lon),length(lat)))
for (i in 1:length(lon)){    
    for (j in 1:length(lat)){
        dat1[i,j] = mean(dat[i,j,])
    }
}
#####################################

library(latticeExtra) # for plotting 
library(maps)         # for ... maps 
library(RColorBrewer)
pal = colorRampPalette(c("blue", "white", "red"))

west.lon.deg=-178.80;
east.lon.deg=178.80; 
step.lon.deg=0.01;
lons  <- c(seq(west.lon.deg, east.lon.deg, step.lon.deg)) 
n.lon <- length(lons) 

south.lat.deg=-90;
north.lat.deg=90;
step.lat.deg=0.01;
lats  <- c(seq(south.lat.deg, north.lat.deg, step.lat.deg)) 
n.lat <- length(lats) 

world.map <- map('world', plot=FALSE, boundary=TRUE,
               xlim=c(west.lon.deg,east.lon.deg), 
               ylim=c(south.lat.deg,north.lat.deg)) 
world.df <- data.frame(lon=world.map$x, lat=world.map$y) 
         
            

l1 = levelplot(c(dat1) ~ coord[,1]*coord[,2], xlab="Longitude",ylab="Latitude", cex.lab=1.2);
l1 + xyplot(lat ~ lon, world.df, type='l', lty=1, lwd=1, col='black', xlab="Longitude",ylab="Latitude", cex.lab=1.2) 

#####################################
medias = c();
desv = c();

for(j in 1:length(lat)){
	medias[j] = mean(dat1[,j]);
	desv[j] = sd(dat1[,j])
}

#####################################
library("mgcv")

ajuste_medias = gam(medias~s(lat,k=20))
ajuste_desv = gam(desv~s(lat,k=20))

lat_new <- seq(-90, 90, length.out = 100)
media_pred <- predict(ajuste_medias, data.frame(lat = lat_new))
desv_pred <- predict(ajuste_desv, data.frame(lat = lat_new))


boxplot(dat1[,order(ncol(dat1):1)],cex=0.2,names=trunc(lat[length(lat):1]),xlab="Latitude",ylab="Temperature",cex.lab=1.2)

plot(lat,medias,xlim=c(-90,90),ylim=c(220,300),xlab="Latitude",ylab="Mean", cex.lab=1.2);
par(new=TRUE)
plot(lat_new,media_pred,type="l",xlim=c(-90,90),ylim=c(220,300),xlab="Latitude",ylab="Mean",lwd = 2, col="mediumorchid3", cex.lab=1.2)

plot(lat,desv,xlim=c(-90,90),ylim=c(-1,20),xlab="Latitude",ylab="Variance", cex.lab=1.2);
par(new=TRUE)
plot(lat_new,desv_pred,type="l",xlim=c(-90,90),ylim=c(-1,20),xlab="Latitude",ylab="Variance",lwd = 2, col="mediumorchid3", cex.lab=1.2)

#####################################
residuos=array(dim=c(length(lon),length(lat)));
for(j in 1:length(lat)){
residuos[,j] = (dat1[,j] - ajuste_medias$fitted.values[j])/(ajuste_desv$fitted.values[j])
}


hist(residuos, col="darkslategray2", main= "Histogram of Residuals", xlab = "Residual", ylab = "Frequency", xlim=c(-10,10), cex.lab=1.2)


l2 = levelplot(residuos ~ coord[,1]*coord[,2], xlab="Longitude",ylab="Latitude", cex.lab=1.2);
l2 + xyplot(lat ~ lon, world.df, type='l', lty=1, lwd=1, col='black', xlab="Longitude",ylab="Latitude", cex.lab=1.2) 

#####################################

x = 6378.388*cos(coord[,2]*pi/180)*cos(coord[,1]*pi/180);
y = 6378.388*cos(coord[,2]*pi/180)*sin(coord[,1]*pi/180);
z = 6378.388*sin(coord[,2]*pi/180);

library("plot3D")
scatter3D(x,y,z,colvar=residuos);


library("gstat");

dataframe = data.frame(cbind(x,y,z,c(residuos)));

g <- gstat(id="residuos", formula=c(residuos)~1, locations=~x+y+z, data = dataframe)
v <- variogram(g,4000,covariogram=TRUE);

plot(v$dist,v$gamma,ylab="Empirical Covariance",xlab="Distance",cex.lab=1.2)
#####################################

ff <- function(alfa){
	suma = 0;
	for(j in 1:length(v$dist)){  
		suma = suma + ( v$gamma[j] - exp(-alfa * (v$dist[j]) ) )^2;    
	}
	return(suma)
	
}

res<-optimize(ff,interval=c(0,0.01));

plot(v$dist,v$gamma,xlim=c(0,4000),ylim=c(0,1.03),ylab="Empirical Covariance",
xlab="Distance",cex.lab=1.2)
par(new=TRUE)
curve(exp(- res$minimum * x),xlim=c(0,4000),ylim=c(0,1.03),col="mediumorchid3",
ylab="Empirical Covariance",xlab="Distance",lwd = 2, cex.lab=1.2)
#####################################

lonmin = -90;
lonmax = -60;
latmin = -60;
latmax = -15;

ind = which(coord[,1]> lonmin & coord[,1]< lonmax & coord[,2]>latmin & coord[,2]<latmax)
N = length(ind);

h = as.matrix(dist(cbind(x[ind],y[ind],z[ind]),diag = TRUE,upper = TRUE));

covmat = array(dim=c(N,N))
for(i in 1:N){
	for(j in 1:N){
		covmat[i,j] = exp(-res$minimum * h[i,j])		
    }
}

inv_covmat = solve(covmat);

#######################################

clon = seq(lonmin,lonmax,l=60);
clat = seq(latmin,latmax,l=60);
coord2 = expand.grid(clon,clat);


cx = 6378.388*cos(coord2[,2]*pi/180)*cos(coord2[,1]*pi/180);
cy = 6378.388*cos(coord2[,2]*pi/180)*sin(coord2[,1]*pi/180);
cz = 6378.388*sin(coord2[,2]*pi/180);

Npred = length(cx);

matpred = array(dim=c(Npred,N));
for(i in 1:Npred){
for(j in 1:N){
	   dij = sqrt( (cx[i] - x[ind][j])^2 + (cy[i] - y[ind][j])^2 + (cz[i] - z[ind][j])^2 );
	   matpred[i,j] = exp(-res$minimum * dij);
}
}	


pred =   matpred%*%inv_covmat%*%cbind(residuos[ind]);


### levelplot datos 
l3 = levelplot(residuos[ind] ~ coord[ind,1]*coord[ind,2],xlab="Longitude",ylab="Latitude", cex.lab=1.2);
l3 + xyplot(lat ~ lon, world.df, type='l', lty=1, lwd=1, col='black', xlab="Longitude",ylab="Latitude", cex.lab=1.2) 


### levelplot downscaling 
l4 = levelplot(c(pred) ~ coord2[,1]*coord2[,2], xlab="Longitude",ylab="Latitude", cex.lab=1.2);
l4 + xyplot(lat ~ lon, world.df, type='l', lty=1, lwd=1, col='black', xlab="Longitude",ylab="Latitude", cex.lab=1.2) 
