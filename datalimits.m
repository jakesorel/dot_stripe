function F = datalimits(mat,clip)

[Nx, Ny, numvar, T] = size(mat);
if T == 1
    mat = cat(4,mat,mat);
end

Mnormvals = zeros(numvar,2);
if clip ~= 0
    for i = 1:numvar
        Mnormvals(i,2) = max(quantile(mat(:,:,i,:),1-clip,'all'));
        Mnormvals(i,1) = min(quantile(mat(:,:,i,:),clip,'all'));
    end
else
    
    for i = 1:numvar
        Mnormvals(i,2) = max(mat(:,:,i,:),[],'all');
        Mnormvals(i,1) = min(mat(:,:,i,:),[],'all');
    end
end
F = Mnormvals;
end
