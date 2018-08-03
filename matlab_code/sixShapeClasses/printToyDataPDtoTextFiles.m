% This is a Matlab script for printing the 6 classes of shape data, generated for example using the command generate_shape_data(n,m,V), to individual text files (on text file for each random sample of points.)

load ToyData_PD_n05.mat
% The cell array ToyData_barcode_n05 is of size 26 x 6 x 2: 25 instances of
% each shape (plus an extra row with shape names at the top), 6 shape
% classes, and 2 homological dimensions (0 and 1).
for i=2:26
    for j=1:6
        for k=1:2
            % We reindix the filenames to get shapes numbered 1-25 and homological dimension numbered 0-1.
            fname = strcat('ToyData_PD_TextFiles/ToyData_PD_n05_',int2str(i-1),'_',int2str(j),'_',int2str(k-1),'.txt');
            fileID = fopen(fname,'w');
            fprintf(fileID,'%.20f %.20f\n',ToyData_barcode_n05{i,j,k}');
            fclose(fileID);
        end
    end
end

load ToyData_PD_n1.mat
for i=2:26
    for j=1:6
        for k=1:2
            % We reindix the filenames to get shapes numbered 1-25 and homological dimension numbered 0-1.
            fname = strcat('ToyData_PD_TextFiles/ToyData_PD_n1_',int2str(i-1),'_',int2str(j),'_',int2str(k-1),'.txt');
            fileID = fopen(fname,'w');
            fprintf(fileID,'%.20f %.20f\n',ToyData_barcode_n1{i,j,k}');
            fclose(fileID);
        end
    end
end
