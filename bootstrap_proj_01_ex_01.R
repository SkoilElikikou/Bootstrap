#Άσκηση 1

data=read.table('D:/OneDrive - unipi.gr/Μεταπτυχιακό/ΥΠΟΛΟΓΙΣΤΙΚΕΣ ΣΤΑΤΙΣΤΙΚΕΣ ΤΕΧΝΙΚΕΣ/ΕΡΓΑΣΙΕΣ/DoctorData.txt',header=TRUE)
attach(data)
n=length(RBCs)
mu_rbc=5
B=1000
alfa=0.05
#α

set.seed(1)
x=replicate(B,RBCs[sample(1:n,n,replace = TRUE)])
xtilde=x-mean(RBCs)
t = (mean(RBCs)-mu_rbc)/(sd(RBCs)/sqrt(n))
t.star = apply(xtilde,2,function(m) mean(m)/(sd(m)/sqrt(n)) )
p=mean(t.star<=t)
p

hist(RBCs,breaks =seq(0,10,0.2))

#b
set.seed(1)
x=replicate(B,RBCs[sample(1:n,n,replace = TRUE)])
x_log=apply(x,c(1,2),log)

x_temp=apply(x_log,2,function(m) sum((m))/n)
x_star=mean(x_temp)
x_star

log_x_cent=x_log-x_star
sigma_temp=apply(log_x_cent,2,function(m) sum(m^2)/n)
sigma_star=mean(sigma_temp)
sigma_star

exp(x_star+ sigma_star/2) #mesh timh ths RBC
(exp(sigma_star)-1)*exp(2*x_star + sigma_star) #diaspora RBC
sqrt((exp(sigma_star)-1)*exp(2*x_star + sigma_star)) #tupikh apoklish RBC


#c

##sunarthsh upologismou tou olikou deikth susxetishs

r_tot=function(data){
  r_temp=NULL
  for (i in 1:(length(data[1,])-1)){
    for (j in (i+1):(length(data[1,]))){
      r_temp=c(r_temp,cor(data[,c(i,j)])[1,2])
    }
  }
  y=mean(r_temp)
  return(y)
}
##########################  San function 
r_total_CIs=function(data,B=1000,alfa=0.05,n=100){
  ##euresh toy gnwstou 'r'
  
  r=r_tot(data)
  data=as.matrix(data)
  set.seed(1)
  y=replicate(B,data[sample(1:n,n,replace = TRUE),])
  
  #euresh tou bootstrap CI
  bootstrap.cint=NULL
  r.star=NULL
  T.star=NULL
  for (b in 1:B){
    r.star[b]=r_tot(y[,,b])
  }
  bootstrap.cint[1]=unname(2*r-rev(quantile(r.star,c(alfa/2,1-alfa/2),type=1)))[1]
  bootstrap.cint[2]=unname(2*r-rev(quantile(r.star,c(alfa/2,1-alfa/2),type=1)))[2]
  bootstrap.cint
  
  #euresh bootstrap-t
  #euresh sd sto kathena bootstrap
  
  bootstrap_t.cint=NULL
  for (b in 1:B){
    
    r.star2=NULL
    for (k in 1:100) {
      l=sample(1:n,n,replace=TRUE)
      r.star2[k]=r_tot(y[,,b][l,])
    }
    T.star[b]=(r.star[b]-r)/sd(r.star2)
  }
  bootstrap_t.cint[1]=r-quantile(T.star,c(1-alfa/2,alfa/2),type=1,names=F)[1]*sd(r.star)
  bootstrap_t.cint[2]=r-quantile(T.star,c(1-alfa/2,alfa/2),type=1,names=F)[2]*sd(r.star)
  
  #euresh quantile bootstrap
  quant_bootstrap.cint=NULL
  quant_bootstrap.cint[1]=quantile(r.star,c(alfa/2,1-alfa/2),type=1,names=F)[1]
  quant_bootstrap.cint[2]=quantile(r.star,c(alfa/2,1-alfa/2),type=1,names=F)[2]
  
  #euresh bca bootstrap
  r.jack = NULL
  bca_bootstrap.cint=NULL
  for ( k in 1:n ) r.jack[k]=r_tot(y[-k,,b])
  a = sum((r.jack-mean(r.jack))^3)/(6*sum((r.jack-mean(r.jack))^2)^(3/2))
  z0 = qnorm(mean(r.star<r))
  z = qnorm(0.975)
  alfa1 = pnorm(z0+(z0-z)/(1-alfa*(z0-z)))
  alfa2 = pnorm(z0+(z0+z)/(1-a*(z0+z)))
  bca_bootstrap.cint[1]=quantile(r.star,c(alfa1,alfa2),type=1,names=F)[1]
  bca_bootstrap.cint[2]=quantile(r.star,c(alfa1,alfa2),type=1,names=F)[2]
  pinakas=matrix(c(bootstrap.cint[1],bootstrap_t.cint[1],quant_bootstrap.cint[1],bca_bootstrap.cint[1],bootstrap.cint[2],bootstrap_t.cint[2],quant_bootstrap.cint[2],bca_bootstrap.cint[2]),ncol = 2)
  row.names(pinakas)=c('Boostrap','Bootstrap-t','Quantile bootstrap','Bca bootstrap')
  colnames(pinakas)=c('Lower','Upper')
  return(pinakas)
}
r_total_CIs(data)



