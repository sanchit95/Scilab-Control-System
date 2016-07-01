//Example 1:
s = %s;
sys = syslin('c',(s+8)/(s^3+4*s^2+5*s+14)) ;
allmargin(sys) 

//Example 2:
z = %z ;
sysd = syslin(0.5,z/(z^3+7*z^2+4*z+8));
a  = allmargin(sys) ;
GainMargin = a.GM ;
PhaseMargin = a.PM ;
