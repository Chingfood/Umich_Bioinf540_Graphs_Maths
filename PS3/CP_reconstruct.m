function temp4 = CP_reconstruct(u1,u2,u3,lambda)
    temp = u1(:,1) * u2(:,1)';
    temp3 = temp * u3(1,1);
    for i = 2:length(u3(:,1))
        temp3 = cat(3,temp3,temp*u3(i,1));
    end
    temp4 = temp3 * lambda(1);
    for j = 2:length(lambda)
        temp = u1(:,j) * u2(:,j)';
        temp3 = temp * u3(1,j);
        for i = 2:length(u3(:,j))
            temp3 = cat(3,temp3,temp*u3(i,j));
        end
        temp3 = temp3 * lambda(j);
        temp4 = temp4 + temp3;
    end
    
    
end