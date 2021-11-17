function x = computenextstate_mytrial4(x,u, paradox,snoise)


x=paradox(8)*x+ paradox(9)*u(:,1)+ paradox(10)*u(:,2)+paradox(11)*u(:,3)+paradox(12)*u(:,4)+snoise;
end
 