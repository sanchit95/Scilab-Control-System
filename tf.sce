//Example 1:
// Continuous time transfer function
sys = tf([5 4 2],[4 2 1 35 6])

//Example 2:
// discrete time transfer function
// Sampling time(Ts = 'd') unspecified
[sys details] = tf([3 8 1],[4 7 5 6 1],'d') 

//Example 3:
[a b]= tf({[3 2 4];[6 4 -9]; [2 3 4]},{[3 4 5];[-3 2 12]; [5 3 4]},0.1) ;
nume = b.num 
deno = b.den
Sampling_Time = b.Ts
name_inp = b.InputName
name_out = b.OutputName



