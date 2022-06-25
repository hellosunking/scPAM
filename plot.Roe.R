#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#
options( stringsAsFactors=F );

outfileName='PAM.revised.Roe.pdf';
pdf( outfileName, width=12, height=6 );
par(mar=c(12,3,0.1,0.1));

Roe = read.table( "Roe.table", head=T );
#added = c(0, 0, 0.5, 1, 1.5, 2.0, 0.1, 1);
#Roe = rbind(Roe, added);

ntype = nrow( Roe );
x = 1:ntype;

enlarge.factor = 3;

##########################################################################################
## Kim
target="Kim";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
plot( x, rep(1,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs,
		ylim=c(0.5, 4.5), axes=F, xlab="", ylab="" );
box();
axis(1, at=1:ntype, cex.lab=1.8, labels=Roe$Celltype, las=2)
axis(2, at=1:4, cex.lab=1.8, labels=c("Kim", "CTR", "PAM", "Legend"), las=2);

## CTR
target="CTR";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
points( x, rep(2,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs);

## PAM
target="PAM";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
points( x, rep(3,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs);

## point-scale
points( 2:12, rep(4,11), type='p', cex=enlarge.factor*c(0.5,1,1.5,0,2, 0,0, 1,0,0,1), pch=16,
	   col=c(rep("grey", 7), "red", "grey", "grey", "blue") );

##########################################################################################
## sc27
target="sc27";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
plot( x, rep(1,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs,
		ylim=c(0.5, 4.5), axes=F, xlab="", ylab="" );
box();
axis(1, at=1:ntype, cex.lab=1.8, labels=Roe$Celltype, las=2)
axis(2, at=1:4, cex.lab=1.8, labels=c("sc27", "sc29", "PAM", "Legend"), las=2);

## sc29
target="sc29";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
points( x, rep(2,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs);

## PAM
target="PAM";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
points( x, rep(3,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs);

## point-scale
points( 2:12, rep(4,11), type='p', cex=enlarge.factor*c(0.5,1,1.5,0,2, 0,0, 1,0,0,1), pch=16,
	   col=c(rep("grey", 7), "red", "grey", "grey", "blue") );

##########################################################################################
## L17
target="L17";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
plot( x, rep(1,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs,
		ylim=c(0.5, 4.5), axes=F, xlab="", ylab="" );
box();
axis(1, at=1:ntype, cex.lab=1.8, labels=Roe$Celltype, las=2)
axis(2, at=1:4, cex.lab=1.8, labels=c("L17", "L20", "CTR", "Legend"), las=2);

## L20
target="L20";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
points( x, rep(2,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs);

## CTR
target="CTR";
cs = c();
for(i in x) {
	if( Roe[i, target] > 1.25 ) {
		cs = c( cs, "red" );
	} else if( Roe[i, target] < 0.8 ) {
		cs = c( cs, "blue" );
	} else {
		cs = c( cs, "grey" );
	}
}
cs = c( cs, "grey" );	## legend
points( x, rep(3,ntype), type='p', cex=enlarge.factor*Roe[, target], pch=16, col=cs);

## point-scale
points( 2:12, rep(4,11), type='p', cex=enlarge.factor*c(0.5,1,1.5,0,2, 0,0, 1,0,0,1), pch=16,
	   col=c(rep("grey", 7), "red", "grey", "grey", "blue") );

