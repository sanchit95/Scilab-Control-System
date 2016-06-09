/////Function :-- stepinfo() ////////////
////Example 1
//Continuous Time System in rational form
s = %s ;
sys = syslin('c',(s+4)/(s^2+8*s+9)) ;
stepinfo(sys) ;

////Example 2:
//Continuous Time System in State-space form
A = [0 1;-5 -2];
B = [0;3];
C = [0 1];
sys_ss = syslin('c',A,B,C) ;
stepinfo(sys_ss,'SettlingTimeThreshold',0.05) ;

////Example 3:
// Discrete Time System
z = poly(0,'z');
sys = syslin('d',z/(z^3+6*z^2+z+6));
stepinfo(sys) ;

////Example 4:
//Assigning values of variables 
y = [2;4;6;8;10;12;14;16;18;20];
t = [1;2;3;4;5;6;7;8;9;10] ;
YFINAL = Y(length(Y));
[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10] = stepinfo(Y,T,YFINAL);
//variables x9 and x10 will give the values of RiseTimeHigh and RiseTimeLow respectively

////Example 5:
//RiseTimeLims
s = %s ;
sys = syslin('c',(s+3)/(2*s^2+4*s+5)) ;
stepinfo(sys,'RiseTimeLims',[0.05 0.95]) ;





