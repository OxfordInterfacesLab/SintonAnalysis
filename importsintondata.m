function DataCell=importsintondata(FileName,PathName,DataCell)

    DataIndex=size(DataCell);
    [~,~,ReadData]=xlsread(strcat(PathName,FileName),'User','A2:O9','basic');
    
    if strcmp(ReadData{1,12}(1),'3') % 3.2 Software version
        num = xlsread(strcat(PathName,FileName),'Calc','I11:L131','basic');
        DataCell(DataIndex(1)+1,1)={[num(:,1),num(:,4)]};
    elseif strcmp(ReadData{1,12}(1),'4') || strcmp(ReadData{1,14}(1),'4')% 4 Software version
        num = xlsread(strcat(PathName,FileName),'RawData','E4:G124','basic');   %% version 4.6 
        DataCell(DataIndex(1)+1,1)={[num(:,3),num(:,1)]};
    end
        %PathName,FileName,SampleName,datenumber, thickness, user resistivity, calc res, 1 sunsIvoc, wafertype
        DataCell(DataIndex(1)+1,2:11)={PathName,FileName,ReadData{5,1},x2mdate(sum(cell2mat(ReadData(5,9:10)))),ReadData{5,2},ReadData{5,3},ReadData{8,3},ReadData{8,11}, ReadData{5,4},...
            mean(DataCell{DataIndex(1)+1,1}((abs(DataCell{DataIndex(1)+1,1}(:,1)-1e15))<1e14,2))};
   
        
end