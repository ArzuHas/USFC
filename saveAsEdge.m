function saveAsEdge(matrix, filename)
    % Save matrix in .edge format with tab delimiters
    fileID = fopen(filename,'w');
    [nRows, nCols] = size(matrix);
    for i = 1:nRows
        for j = 1:nCols
            if j == nCols
                fprintf(fileID, '%f', matrix(i,j));
            else
                fprintf(fileID, '%f\t', matrix(i,j));  % Use tab as delimiter
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end
