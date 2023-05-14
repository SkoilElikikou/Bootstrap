data=read.table('D:/OneDrive - unipi.gr/Μεταπτυχιακό/ΥΠΟΛΟΓΙΣΤΙΚΕΣ ΣΤΑΤΙΣΤΙΚΕΣ ΤΕΧΝΙΚΕΣ/ΕΡΓΑΣΙΕΣ/HeightsData.txt',header=TRUE)
attach(data)
n=length(Height)
euros=function(x) {max(x)-min(x)}

jknife=function(x,fxn,alfa=0.05) {
  fake_x=as.matrix(x)
  if (ncol(as.matrix(x))=='ΝΑ') {
    
  theta=fxn(x)
  n=length(x)
  partials=NULL
  for (i in 1:n) {
    partials[i]=fxn(x[-i])
  }
  pseudos=(n*theta)-(n-1)*partials
  jack.bias=mean((n-1)*(partials-theta))
  jack.est=mean(pseudos)
  jack.se=sd(pseudos)/sqrt(n)
  CI=qt(alfa/2,n-1,lower.tail = FALSE)*jack.se
  jack.ci=c(jack.est-CI,jack.est+CI)
  list(est=jack.est,se=jack.se,ci=jack.ci,pseudos=pseudos,partials=partials,bias=jack.bias)
  }
  else {
    theta=fxn(x)
    n=nrow(x)
    partials=NULL
    for (i in 1:n) {
      partials[i]=fxn(x[-i,])
    }
    pseudos=(n*theta)-(n-1)*partials
    jack.bias=mean((n-1)*(partials-theta))
    jack.est=mean(pseudos)
    jack.se=sd(pseudos)/sqrt(n)
    CI=qt(alfa/2,n-1,lower.tail = FALSE)*jack.se
    jack.ci=c(jack.est-CI,jack.est+CI)
    list(est=jack.est,se=jack.se,ci=jack.ci,pseudos=pseudos,partials=partials,bias=jack.bias)
  }
  }
jknife(Height,euros) 

n*jknife(Height,euros)$theta-(n-1)*jknife(Height,euros)$asdf

S=0
for (i in 1:n) {
  if (data$Height[i]>165) S=S+1
}
S/n
mean(data$Height>165)

jknife(Height,sd)
sd(Height)


####Δ.Ε.
bootstrap_CIs(Height,fxn=sd,fxn_dim = 1,alfa = 0.02)
sd(Height)
sd(Height)+c(-1,1)*qnorm(0.98)*0.2577
