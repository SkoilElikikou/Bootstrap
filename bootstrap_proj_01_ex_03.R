data=c(12,24,16,28,17,19,31,15,38,19,8,20,12,21,14,15,22,20,24,14,21,25,18,27,18,16,29,14,24,14,10,31,11,23,24,8,29,10,23,23,19,29,13,21,23,16,37,10,29,23)
matrix(data,nrow=5)
n=length(data)


####mesh diafora####
win_diff=function(data) {
  n=length(data)
  desc_data=sort(data,decreasing = TRUE)
  s_first=sum(desc_data[1:10])
  s_last=sum(desc_data[(n-9):n])
  return(s_first-s_last)
}
win_diff(data)

#####diafora meswn#####
diff_win=function(data) {
  n=length(data)
  desc_data=sort(data,decreasing = TRUE)
  s_first=desc_data[1:10]
  s_last=desc_data[(n-9):n]
  return(mean(s_first)-mean(s_last))
}
diff_win(data)
ncol(as.matrix(data))
bootstrap_CIs(data,fxn=diff_win,alfa=0.05,fxn_dim = 1)
bootstrap_CIs(data,fxn=win_diff,alfa=0.05,fxn_dim = 1)

sum(data[1:10])-sum(data[41:50])
data=as.vector(data)

sort(data,decreasing = TRUE)


hist(data,breaks=seq(0,50,1))




#verifications
set.seed(1)
y=replicate(1000,data[sample(1:n,n,replace=TRUE)])
y
dim(y)
y[,1]
for (i in 1:1000) {
  cat(win_diff(y[,i]),'\n')
}
win_diff(y[,2])
win_diff(y[,3])
mean(win_diff(y))
y_star=apply(y,2,win_diff)
mean(y_star)
sd(y_star)
mean(y_star)+c(-1,1)*qnorm(0.95)*sd(y_star)
