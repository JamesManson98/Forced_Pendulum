function Gu=Gu(y,s,j,M)
df1=@(u,a)MyJacobian(@(u)M(u,a),u,1e-5);    %we define df1, differentiating M with respect to u
Gu=cell(2^(j-1),2^(j-1));   %we preallocate the cell Gu
for x=1:2^(j-1)
    for z=1:2^(j-1)
        Gu{x,z}=zeros(2);   %we fill the cell with 2x2 arrays of zeros
    end
end
for k=1:2^(j-1)
    J=df1(y(2*k-1:2*k),y(end)); %we find df1 of a point in y (partial M/partial u)
    Gu{k,k}=J;  %we set the diagonal as J1,J2,J3..Jj
    if k<2^(j-1)
        Gu{k,k+1}=(-1)*eye(2);  %we set the matrix one to the right of J1,...,Jj as the negative of the 2x2 identity matrix
    end
end
if j>1
    Gu{2^(j-1),1}=s*eye(2); %the bottom left matrix is defined as s*the 2x2 identity matrix
end
Gu=cell2mat(Gu);    %we then convert this back into a matrix to output
end