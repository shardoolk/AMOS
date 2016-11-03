function v = f1(l,l1,l2)
global K 
global Q 
v = 1.5*(l2*l2)/l1 - (l1/(1-Q*l))*(2*K - Q*l2);
