function Ga=Ga(y,j,M)
df1=@(u,a)MyJacobian(@(a)M(u,a),a,1e-6);    %we define df1 as partial M by partial a
Ga=cell(2^(j-1),1); %we preallocate the cell Ga
for i=1:2^(j-1)
    Ga{i,1}=df1(y(2*i-1:2*i),y(end));   %we put the pairs of y values into df1 and put them into Ga
end
Ga=cell2mat(Ga);    %we convert this back into a matrix, Ga, to output.
end