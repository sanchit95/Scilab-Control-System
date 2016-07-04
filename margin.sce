//Example 1:
s = %s ;
sys = syslin('c',3/(s^3+7*s^2+8*s+5)) ;
[GM,PM,Fgm,Fpm] = margin(sys) ;

//Example 2:
z = %z ;
sys1 = syslin('d',4/(z^2+7*z+8)) ;
margin(sys1) ;  //Plots bode plot showing phase margin, gain margin etc. 
