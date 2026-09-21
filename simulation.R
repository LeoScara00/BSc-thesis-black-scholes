# Reproducible Black-Scholes Monte Carlo study. Base R only.
# Run: Rscript simulation.R [output-directory]
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else 'results'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'figures'), showWarnings=FALSE)
RNGkind('Mersenne-Twister','Inversion','Rejection')
r <- .0329; sigma <- .25; T <- 90/365; s0 <- 24; K <- 26; N <- 400L
cdf <- pnorm
d1 <- (log(s0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T)); d2 <- d1-sigma*sqrt(T)
benchmark <- s0*cdf(d1)-K*exp(-r*T)*cdf(d2)
pay <- function(s) exp(-r*T)*pmax(s-K,0)
summary_mean <- function(y) {
 est <- mean(y); se <- sd(y)/sqrt(length(y))
 c(estimate=est,se=se,lower=est-1.96*se,upper=est+1.96*se)
}
# Main experiment: independent paths, paired methods, nested 500/1000 samples.
set.seed(20260921)
dw <- matrix(rnorm(1000*N,sd=sqrt(T/N)),nrow=1000,ncol=N)
exact <- em <- matrix(s0,1000,N+1)
for(j in seq_len(N)) {
 exact[,j+1] <- exact[,j]*exp((r-sigma^2/2)*T/N+sigma*dw[,j])
 em[,j+1] <- em[,j]*(1+r*T/N+sigma*dw[,j])
}
stopifnot(max(abs(exact[,N+1]-s0*exp((r-sigma^2/2)*T+sigma*rowSums(dw))))<1e-9)
main <- do.call(rbind,lapply(c(500,1000),function(m) do.call(rbind,lapply(c('Exact','Euler'),function(method) {
 y <- pay(if(method=='Exact') exact[seq_len(m),N+1] else em[seq_len(m),N+1])
 data.frame(method=method,paths=m,as.list(summary_mean(y)),absolute_error=abs(mean(y)-benchmark))
}))))
write.csv(main,file.path(out,'main_results.csv'),row.names=FALSE)
write.csv(data.frame(path=1:1000,exact_terminal=exact[,N+1],euler_terminal=em[,N+1],exact_payoff=pay(exact[,N+1]),euler_payoff=pay(em[,N+1])),file.path(out,'main_samples.csv'),row.names=FALSE)
# Exact sampling: independent replications; nested sample sizes within each replication.
set.seed(20260922)
sizes <- c(100,500,1000,5000,10000); B <- 200L
estimates <- matrix(NA_real_,B,length(sizes)); covers <- estimates
for(b in seq_len(B)) {
 y <- pay(s0*exp((r-sigma^2/2)*T+sigma*sqrt(T)*rnorm(max(sizes))))
 for(j in seq_along(sizes)) {
  q <- summary_mean(y[seq_len(sizes[j])]); estimates[b,j] <- q[1]
  covers[b,j] <- q[3]<=benchmark && benchmark<=q[4]
 }
}
mc <- data.frame(paths=sizes,replications=B,mean_price=colMeans(estimates),bias=colMeans(estimates)-benchmark,rmse=sqrt(colMeans((estimates-benchmark)^2)),empirical_sd=apply(estimates,2,sd),coverage=colMeans(covers))
# Analytical variance of discounted call payoff, from truncated lognormal moments.
a <- log(s0)+(r-sigma^2/2)*T; v <- sigma^2*T
moment <- function(p) exp(p*a+.5*p*p*v)*pnorm((a+p*v-log(K))/sqrt(v))
varY <- exp(-2*r*T)*(moment(2)-2*K*moment(1)+K*K*moment(0))-benchmark^2
mc$theoretical_se <- sqrt(varY/sizes)
write.csv(mc,file.path(out,'sampling_convergence.csv'),row.names=FALSE)
write.csv(data.frame(replication=rep(seq_len(B),times=length(sizes)),paths=rep(sizes,each=B),price=as.vector(estimates)),file.path(out,'sampling_replications.csv'),row.names=FALSE)
# Coupled grid study: aggregate fine increments; same Brownian paths at every resolution.
# Pairwise differences isolate discretisation from the much larger payoff variance.
set.seed(20260923)
M <- 20000L; steps <- c(25L,50L,100L,200L,400L)
fine <- matrix(rnorm(M*N,sd=sqrt(T/N)),M,N)
terminal <- s0*exp((r-sigma^2/2)*T+sigma*rowSums(fine))
differences <- matrix(NA_real_,M,length(steps))
grid <- do.call(rbind,lapply(seq_along(steps),function(k) {
 n <- steps[k]; R <- N/n; s <- rep(s0,M); ever_negative <- rep(FALSE,M)
 for(j in seq_len(n)) {
  increment <- rowSums(fine[,((j-1)*R+1):(j*R),drop=FALSE])
  s <- s*(1+r*T/n+sigma*increment); ever_negative <- ever_negative | s<0
 }
 differences[,k] <<- pay(s)-pay(terminal)
 data.frame(steps=n,dt=T/n,paths=M,as.list(summary_mean(differences[,k])),terminal_rmse=sqrt(mean((s-terminal)^2)),negative_path_fraction=mean(ever_negative))
}))
write.csv(grid,file.path(out,'euler_convergence.csv'),row.names=FALSE)
write.csv(data.frame(path=seq_len(M),setNames(as.data.frame(differences),paste0('difference_N',steps))),file.path(out,'paired_payoff_differences.csv'),row.names=FALSE)
# Both PDF (LaTeX) and PNG (preview). No external plotting packages.
figure <- function(name,draw,width=8,height=5) {
 pdf(file.path(out,'figures',paste0(name,'.pdf')),width=width,height=height,useDingbats=FALSE); draw(); dev.off()
 png(file.path(out,'figures',paste0(name,'.png')),width=width*160,height=height*160,res=160); draw(); dev.off()
}
figure('trajectories',function() {
 par(mfrow=c(1,2),mar=c(4,4,3,1))
 for(method in c('Exact','Euler')) {
  x <- if(method=='Exact') exact else em
  matplot(seq(0,T,length.out=N+1),t(x[1:20,]),type='l',lty=1,col=adjustcolor('#176b91',alpha.f=.55),xlab='Time (years)',ylab='Stock price',main=paste(method,': 20 paired paths'),ylim=range(exact[1:20,],em[1:20,]))
 }
},10,4.5)
figure('terminal_distributions',function() {
 par(mfrow=c(1,2),mar=c(4,4,3,1)); breaks <- pretty(range(exact[,N+1],em[,N+1]),20)
 for(method in c('Exact','Euler')) {
  x <- if(method=='Exact') exact[,N+1] else em[,N+1]
  hist(x,breaks=breaks,probability=TRUE,col='#c6dce7',border='white',main=paste(method,': 1,000 paths'),xlab='Terminal stock price',ylim=c(0,.16))
  curve(dlnorm(x,meanlog=a,sdlog=sqrt(v)),add=TRUE,col='#ab432e',lwd=2)
  legend('topright','Exact lognormal density',col='#ab432e',lty=1,bty='n',cex=.75)
 }
},10,4.5)
figure('sampling_convergence',function() {
 par(mfrow=c(1,2),mar=c(4,4,3,1))
 plot(mc$paths,mc$rmse,log='xy',type='b',pch=19,col='#176b91',xlab='Number of paths',ylab='Pricing RMSE',main='200 independent replications')
 lines(mc$paths,mc$theoretical_se,col='#ab432e',lwd=2,lty=2)
 legend('topright',c('Empirical RMSE','Theoretical standard error'),col=c('#176b91','#ab432e'),lty=c(1,2),bty='n',cex=.75)
 boxplot(as.data.frame(estimates),names=sizes,xlab='Number of paths',ylab='Estimated call price',main='Distribution across replications',col='#c6dce7',outline=FALSE)
 abline(h=benchmark,col='#ab432e',lwd=2,lty=2)
},10,4.5)
figure('euler_convergence',function() {
 par(mfrow=c(1,2),mar=c(4,4.5,3,1))
 plot(grid$dt,grid$terminal_rmse,log='xy',type='b',pch=19,col='#176b91',xlab='Time step (years)',ylab='Terminal stock RMSE',main='Paired path approximation')
 lines(grid$dt,tail(grid$terminal_rmse,1)*sqrt(grid$dt/tail(grid$dt,1)),col='#ab432e',lty=2)
 legend('topleft',c('Observed RMSE','Square-root reference'),col=c('#176b91','#ab432e'),lty=c(1,2),bty='n',cex=.75)
 plot(grid$steps,grid$estimate,log='x',pch=19,col='#176b91',ylim=range(grid$lower,grid$upper,0),xlab='Number of time steps',ylab='Euler minus exact call price',main='Paired mean and 95% interval')
 arrows(grid$steps,grid$lower,grid$steps,grid$upper,angle=90,code=3,length=.05,col='#176b91');abline(h=0,lty=2,col='#ab432e')
},10,4.5)
stopifnot(all(is.finite(main$estimate)),all(is.finite(grid$terminal_rmse)),length(unique(em[,N+1]))==1000)
writeLines(c(sprintf('Black-Scholes benchmark: %.12f',benchmark),sprintf('Analytical discounted-payoff variance: %.12f',varY),'Seeds: 20260921 (main), 20260922 (sampling), 20260923 (Euler).','Intervals are asymptotic marginal 95% normal intervals, not simultaneous intervals.','Euler grid estimates share paths; sampling sizes are nested within each independent replication.',capture.output(sessionInfo())),file.path(out,'run_metadata.txt'))
print(main);print(mc);print(grid)
# LaTeX tables generated from the same result objects (no manual transcription).
main_lines <- vapply(seq_len(nrow(main)),function(i) with(main[i,],sprintf('%s & %d & %.5f & %.5f & [%.5f, %.5f] \\\\',method,paths,estimate,se,lower,upper)),character(1))
writeLines(c('\\begin{tabular}{lrrrr}','\\toprule','Method & Paths & Price & Std. error & 95\\% interval\\\\','\\midrule',main_lines,'\\bottomrule','\\end{tabular}'),file.path(out,'main_table.tex'))
grid_lines <- vapply(seq_len(nrow(grid)),function(i) with(grid[i,],sprintf('%d & %.5f & %.6f & [%.6f, %.6f] \\\\',steps,terminal_rmse,estimate,lower,upper)),character(1))
writeLines(c('\\begin{tabular}{rrrr}','\\toprule','Steps & Stock RMSE & Price difference & 95\\% interval\\\\','\\midrule',grid_lines,'\\bottomrule','\\end{tabular}'),file.path(out,'euler_table.tex'))
