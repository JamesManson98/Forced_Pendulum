function [ynew,xnew,snew]=MapLine(M,x,s,maxdist,maxangle)
xnew=x; %we find the initial xnew
snew=s; %we find the initial snew
ynew=M(xnew);   %we find ynew by mapping all of the points in xnew to ynew
i=2;    
while i<=length(ynew)   %we go through every point in ynew
    while norm(ynew(:,i)-ynew(:,i-1))>maxdist   %we find where the distance between points is too great
        s_star=(snew(i)+snew(i-1))/2;   %we find a new point to insert into x
        x_star = interp1(snew,xnew',s_star,'spline');   %we find a new x using interp1
        snew=[snew(1:i-1),s_star,snew(i:end)];  %we insert the new point into snew
        xnew=[xnew(:,1:i-1),x_star',xnew(:,i:end)]; %we insert the new point into xnew
        ynew=[ynew(:,(1:i-1)),M(x_star'),ynew(:,(i:end))];  %we calculate and insert the equivalent point for y into ynew
    end
    i=i+1;  %we add 1 to i and check whether there needs to be a point added between the new i and i-1
end   
i=2;
while i<length(ynew)    %we go through every point in ynew
    while abs(abs((180/pi)*acos((dot(ynew(:,i+1)-ynew(:,i),ynew(:,i-1)-ynew(:,i))/(norm(ynew(:,i+1)-ynew(:,i))*norm(ynew(:,i-1)-ynew(:,i))))))-180)>maxangle    %we find where the angles are too large and a point must be inserted
            s_star=(snew(i)+snew(i-1))/2;   %we find the 2 points to insert into snew
            s_star1=(snew(i)+snew(i+1))/2;
            x_star=interp1(snew,xnew',s_star,'spline');     %we use interp1 to find the 2 points to insert into x_new
            x_star1=interp1(snew,xnew',s_star1,'spline');
            snew=[snew(1:i-1),s_star,snew(i),s_star1,snew(i+1:end)];                    %we insert the corresponding points into ynew, xnew, and snew
            xnew=[xnew(:,1:i-1),x_star',xnew(:,i),x_star1',xnew(:,i+1:end)];
            ynew=[ynew(:,(1:i-1)),M(x_star'),ynew(:,i),M(x_star1'),ynew(:,(i+1:end))];
    end
    i=i+1;  %we add 1 to i and check whether there needs to be a point added between the new i and i-1
end
i=3;
while i<length(ynew)-1
    if norm(ynew(:,i+1)-ynew(:,i-1))<maxdist    %we check whether we can remove a point if the distance between i+1 and i-1 are satisfied
    if abs(abs((180/pi)*acos((dot(ynew(:,i+1)-ynew(:,i-1),ynew(:,i-2)-ynew(:,i-1))/(norm(ynew(:,i+1)-ynew(:,i-1))*norm(ynew(:,i-2)-ynew(:,i-1))))))-180)<maxangle&&abs(abs((180/pi)*acos((dot(ynew(:,i+2)-ynew(:,i+1),ynew(:,i-1)-ynew(:,i+1))/(norm(ynew(:,i+2)-ynew(:,i+1))*norm(ynew(:,i-1)-ynew(:,i+1))))))-180)<maxangle   %we check whether both angles that involve i will be satisfied with i removed
            snew=[snew(1:i-1),snew(i+1:end)];       %we delete the point snew(i) from snew if both conditions are satisfied
            xnew=[xnew(:,1:i-1),xnew(:,i+1:end)];   %we delete the point xnew(:,i) from xnew if both conditions are satisfied
            ynew=[ynew(:,1:i-1),ynew(:,i+1:end)];   %we delete the point ynew(:,i) from ynew if both conditions are satisfied
            continue
    end
    end
    i=i+1;
end
if norm(ynew(:,3)-ynew(:,1))<maxdist %we test whether we can remove the point ynew(:,2) by seeing whether the conditions are satisfied between points 1, 3 and 4
    if abs(abs((180/pi)*acos((dot(ynew(:,4)-ynew(:,3),ynew(:,1)-ynew(:,3))/(norm(ynew(:,4)-ynew(:,3))*norm(ynew(:,1)-ynew(:,3))))))-180)<maxangle
            snew=[snew(1),snew(3:end)];
            xnew=[xnew(:,1),xnew(:,3:end)];
            ynew=[ynew(:,1),ynew(:,3:end)]; 
    end
end  
if norm(ynew(:,end)-ynew(:,end-2))<maxdist %we test whether we can remove the point ynew(:,end-1) by seeing whether the conditions are satisfied between points end-3, end-2 and end.
    if abs(abs((180/pi)*acos((dot(ynew(:,end)-ynew(:,end-2),ynew(:,end-3)-ynew(:,end-2))/(norm(ynew(:,end)-ynew(:,end-2))*norm(ynew(:,end-3)-ynew(:,end-2))))))-180)<maxangle
            snew=[snew(1:end-2),snew(end)];
            xnew=[xnew(:,1:end-2),xnew(:,end)];
            ynew=[ynew(:,1:end-2),ynew(:,end)]; 
    end
end
end