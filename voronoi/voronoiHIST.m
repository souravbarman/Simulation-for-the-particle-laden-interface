function voronoiHIST
    matfiles=dir(fullfile('*.tif'));
    totalFiles=length(matfiles);
    disp(totalFiles);
    finalTally=zeros(20,2);
    for k=1:20
        finalTally(k,1)=k;
    end
    for i=1:totalFiles
        disp(matfiles(i).name);
%         test=i/50;
%         if isinteger(test)
%             disp(i);
%         end
        currentImage=imread(matfiles(i).name);
        [row,col]=find(currentImage==0);
        points=[col,row];

        try
            [v,c]=voronoin(points);
        catch err
            continue
        end
        
        totalCells=length(c);
        symetry=zeros(totalCells,1);
        for j=1:totalCells
            symetry(j,1)=length(c{j});
        end
        for k=1:20
            finalTally(k,2)=finalTally(k,2)+length(find(symetry==k));
        end
        
    end
    xlswrite('voronoiHIST_results.xls', finalTally);
end
    