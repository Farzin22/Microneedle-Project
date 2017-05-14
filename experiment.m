%experiment analysis matlab code
%author Farzin Akbar


number_of_experiments=1  %Only one experiment at a time is analysed

%g is for figure(2)
g = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(g, 'Position'); % grab original on screen size
set(g, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(g,'PaperPositionMode','auto') %set paper pos for printing



for j=1:number_of_experiments 
    
    
    
    %Estimate the surface meeting force 
    val_forcy =force_threshold; %value to find
    forcies= tsmovavg(FDs([1:2048],2),'t',FD_movingaverage,1);
    tmp1 = abs(forcies-val_forcy)
    [idxx idxx] = min(tmp1) %index of closest value
    closest_force = forcies(idxx) %closest value
    index=idxx; %finding the index of the threshold voltage
    initial_distance_mm =  FDs(index,1)-FDs(1,1); 
	initial_time_s(j) = FTs(index,1)-FTs(1,1);
    initial_distance(j) = initial_distance_mm * 1000; %µm
    
	%calculating the differentiate of the Force-Index diagram in order to find the penetration force
	begin_diff=100;
	end_diff=2048;
    subplot(2,2,1)
    displacement_offset(j) = FDs(1,1);
    dydx = diff([eps;tsmovavg(FDs([1:2048],2),'t',FD_movingaverage,1)*1000])./diff([eps;(FDs([1:2048],1)-displacement_offset(j))*1000-initial_distance(j)]);
    xc=1:1:2048;
    [M,I] = min(dydx(begin_diff:end_diff));
	
	while(abs(M)>10 )
		
		end_diff=end_diff-100;
		if(end_diff < begin_diff)
			end_diff=end_diff+100;
		end
		[M,I] = min(dydx(begin_diff:end_diff));
	end
    I=I+begin_diff-1; %to find the real index!
    b(j)= plot((FDs([1:2048],1)-displacement_offset(j))*1000-initial_distance(j),tsmovavg(FDs([1:2048],2),'t',FD_movingaverage,1)*1000,'-','Linewidth',2);
    hold on
	popo=tsmovavg(FDs(:,2),'t',FD_movingaverage,1)*1000;
    
    %calculation of the penetration force, displacementat penetration, drope of the force at penetration, time of the penetration 
    penetration_force(j)=max(popo(1:I));
    dropped_force (j) = min(popo(I:length(popo)));
    penetration_force_index(j)= find(popo==penetration_force(j));
    dropped_force_index(j) = min(find(popo==dropped_force(j)));
    force_drop(j)= penetration_force(j)- dropped_force (j);
	%plotting the force-displacement, force-time, velocity-displacement, and displacement-time diagrams.
    displacement_at_penetration(j)= (FDs(penetration_force_index(j),1)-displacement_offset(j))*1000-initial_distance(j) ;
    plot((FDs(penetration_force_index(j),1)-displacement_offset(j))*1000-initial_distance(j),popo(penetration_force_index(j)),'*');
    plot((FDs(dropped_force_index(j),1)-displacement_offset(j))*1000-initial_distance(j),popo(dropped_force_index(j)),'*');
    set(gca,'fontsize',20)
    xlabel('Displacement (µm)');
    ylabel('Force (mN)');
    grid on

    hold on

    subplot(2,2,2);
    
    
    c(j)=plot(FTs([1:2048],1)-initial_time_s(j),(tsmovavg(FTs([1:2048],2),'t',FT_movingaverage,1))*1000,'-','Linewidth',2);
    time_at_penetration(j)= FTs(penetration_force_index(j),1)-initial_time_s(j);
    set(gca,'fontsize',20)
    
    %ylim([0  0.2])
    hold on

    xlabel('Time (s)');
    ylabel('Force (mN)');
    
    hold on
    grid on
    
    subplot(2,2,3);
    Vdisplacement_offset(j) = VDs(1,1);
    d(j)= plot((VDs([1:2048],1)-Vdisplacement_offset(j))*1000-initial_distance(j),tsmovavg(VDs([1:2048],2),'t',VD_movingaverage,1),'-','Linewidth',2);
    set(gca,'fontsize',20)
    xlabel('Displacement (µm)');
    ylabel('Velocity (mm/s)');
    hold on
    grid on
    
    subplot(2,2,4);
    Ddisplacement_offset(j) = DTs(1,j);
    e(j)=plot(FTs([1:2048],1)-initial_time_s(j),(tsmovavg(DTs([1:2048],2),'t',DT_movingaverage,1)-Ddisplacement_offset(j))*1000-initial_distance(j),'-','Linewidth',2);
    set(gca,'fontsize',20)
    xlabel('Time (s)');
    ylabel('Displacement (µm)');
    hold on
    grid on
    
end
%writing titiles of the plots
subplot(2,2,1)
str50 = ['$N_\chi$=', sprintf('%d, ', FD_movingaverage)];% , '$\sigma$=', ...
str60 = ['$F_{penetration}$ =', sprintf('%5.2f mN, ', mean(penetration_force)), '$F_{drop}$ =', sprintf('%5.2f mN, ', mean(force_drop))];
str70 = ['$D_{penetration}$ =', sprintf('%5.2f ', mean(displacement_at_penetration)),'$\mu$ ' ,'$m$'] ;
str_final=strvcat(str50,str60,str70);
title(str_final,'interpreter','latex');
grid on
grid minor

subplot(2,2,2);
str500 = ['$N_\chi$=', sprintf('%d, ', FT_movingaverage)];
str600 = ['$T_{penetration}$ =', sprintf('%5.3f s, ', mean(time_at_penetration))];
str_final_time=strvcat(str500,str600);
title(str_final_time,'interpreter','latex');
grid on
grid minor

subplot(2,2,3);
grid on
grid minor
str_VD = ['$N_\chi$=', sprintf('%d, ', VD_movingaverage)];
title(str_VD,'interpreter','latex');

subplot(2,2,4);
grid on
grid minor
str_DT = ['$N_\chi$=', sprintf('%d, ', DT_movingaverage)];
title(str_VD,'interpreter','latex');


%%%%%%%%%%%%%%%%%%%%%%%saving the figure %%%%%%%%%%%%%
cd(path);
file_name =strcat(path,'\','analyzed');
print(g,'analyzed','-dpng', '-r300')
print(g,'analyzed','-dtiff')
saveas(g, 'analyzed') % save figure
set(g,'Position', origSize) %set back to original dimensions

%openfig('analyzed.fig','visible')  %write this code of line in matlab if you want to open the figure
%%%%%%%%%reading experiment details from readme.txt%%%%%%%%%%%
read_filename ='readme.txt';
the_path =strcat(read_filename);
dd = dir(read_filename);                           
fileName = {dd.name}; 
file    = memmapfile( char(fileName), 'writable', true );
comma   = uint8(',');
point   = uint8('.');
file.Data( transpose( file.Data==comma) ) = point;
[parameters, values] = textread(read_filename,'%s %s', 11);
IndexC = strfind(parameters, 'Diameter');
Index = find(not(cellfun('isempty', IndexC)));
diameter = str2double(values(Index));
IndexC = strfind(parameters, 'Number');
Index = find(not(cellfun('isempty', IndexC)));
number_of_microneedles = str2double(values(Index));
IndexC = strfind(parameters, 'mm/s');
Index = find(not(cellfun('isempty', IndexC)));
Velocity = str2double(values(Index));
to_print=strcat(file_name,'.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%writing the experiment results into file%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir('Experiment_Analysis')
new_path =strcat(path,'\','Experiment_Analysis');
cd(new_path);

final_analysis_filename= strcat('final_analysis_',num2str(number_of_microneedles),'mn_','V_',num2str(Velocity),'mms-1','.csv');
a=exist(final_analysis_filename, 'file')

if(a==0)
	fid = fopen(final_analysis_filename, 'a');
	fprintf(fid, 'Number of Microneedles ; Tip Diameter (µm); Penetration Force (mN); Force Drop (mN); Displacement at Penetration(µm);Time at Penetration (s);Threshold Force (mN);FD_movingaverage;FT_movingaverage;DT_movingaverage;VD_movingaverage');
	
		
	switch char(parameters(1))
		case 'Position_Controlled'
			fprintf(fid,'; Applied Position (µm)')  
		case 'Velocity_Controlled'
			fprintf(fid,'; Applied Velocity (mm/s)') 
		case 'Force_Controlled'
			fprintf(fid,'; Applied Force (mN)') 
		otherwise
			disp('incorrect mode')
	end
	
	fclose(fid)
end

fid = fopen(final_analysis_filename, 'a');

fprintf(fid,strcat('\n',strrep(num2str(number_of_microneedles),'.', ','),';',strrep(num2str(diameter),'.', ','),';',strrep(num2str(mean(penetration_force)),'.', ',') ,';', strrep(num2str(mean(force_drop)) ,'.', ','),';' ,  strrep(num2str(mean(displacement_at_penetration)),'.', ','),';', strrep(num2str(mean(time_at_penetration)) ,'.', ','),';',strrep(num2str(mean(force_threshold*1000)) ,'.', ','),';',strrep(num2str(mean(FD_movingaverage)) ,'.', ','),';',strrep(num2str(mean(FT_movingaverage)) ,'.', ','),';',strrep(num2str(mean(DT_movingaverage)) ,'.', ','),';',strrep(num2str(mean(VD_movingaverage)) ,'.', ',')));
switch char(parameters(1))
    case 'Position_Controlled'
        fprintf(fid,strcat(';', strrep(num2str(applied_variable*1000),'.', ',') ));  
    case 'Velocity_Controlled'
        fprintf(fid,strcat(';', strrep(num2str(applied_variable),'.', ',') )); 
    case 'Force_Controlled'
        fprintf(fid,strcat(';', strrep(num2str(applied_variable*1000),'.', ',') )); 
    otherwise
        disp('incorrect mode')
end

fclose(fid);


