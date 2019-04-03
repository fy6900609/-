% zrcadleni, Perturbation y into <a,b>
function result=zrcad_shade(y,xi,a,b)
delka=length(y);
for i=1:delka
    if (y(i)<a)||(y(i)>b)
		if y(i)>b
		    y(i)=(b+xi(1,i))/2;
		elseif y(i)<a
		    y(i)=(a+xi(1,i))/2;
		end
    end
end
result=y;