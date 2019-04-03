% zrcadleni, Perturbation y into <a,b>
function result=zrcad(y,a,b)
zrc=find(y<a|y>b);
for i=zrc
	while (y(i)<a|y(i)>b)
		if y(i)>b
		    y(i)=2*b-y(i);
		elseif y(i)<a
		    y(i)=2*a-y(i);
		end
	end
end
result=y;