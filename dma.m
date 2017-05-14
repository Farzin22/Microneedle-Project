%DMA Analysis Code
%Autor Farzin Akbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DMA probe:
probe_diameter_micrometer = 8346;
probe_radius = (probe_diameter_micrometer/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
my_sweep_index=num2str(sweep_index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(f, 'Position'); % grab original on screen size
set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(f,'PaperPositionMode','auto') %set paper pos for printing

%f is for figure(1)    
subplot(2,2,1)
b = plot(FDs([1:2048],1)*1000,tsmovavg(FDs([1:2048],2),'t',FD_movingaverage,1),'-');
set(gca,'fontsize',20)
xlabel('Displacement (µm)','FontSize', 20);
ylabel('Force (N)','FontSize', 20);
grid on
grid minor    
hold on
    
subplot(2,2,2);
c=plot(FTs([1:2048],1),tsmovavg(FTs([1:2048],2),'t',FT_movingaverage,1),'-');
set(gca,'fontsize',20)
xlabel('Time (s)','FontSize', 20);
ylabel('Force (N)','FontSize', 20);
grid on
grid minor
hold on

%subplot(2,2,3);
%hold on
%grid on

subplot(2,2,4);
plot(FTs([1:2048],1),tsmovavg(DTs([1:2048],2),'t',DT_movingaverage,1)*1000,'-');
set(gca,'fontsize',20)
xlabel('Time (s)','FontSize', 20);
ylabel('Displacement (µm)','FontSize', 20); 
hold on
grid minor
grid on

file_name =strcat(path,'\',Ex_name,'_',my_sweep_index,'_raw_data')

print(f,file_name,'-dpng', '-r300')
print(f,file_name,'-dtiff')
saveas(f, file_name) % save figure



set(f,'Position', origSize) %set back to original dimensions

to_print_raw=strcat(file_name,'.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%Estimate the steady state time 
timies= FTs([1:2048],1);
tmp1 = abs(timies-val_timy)
[idxx idxx] = min(tmp1) %index of closest value
closest_time = timies(idxx) %closest value
index=idxx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

indentation_Ds = DTs ([index:2048],2);
indentation_Fs = FTs ([index:2048],2);
indentation_ts = FTs ([index:2048],1);

for i =1:2048
   Indentation_depth(i) = (DTs(i,2)-DTs(1,2))*1000;
end

%Finding Zero Crossings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
displacement_peak_to_peak = max (indentation_Ds)- min (indentation_Ds);
first_displacement= indentation_Ds(1);

displaced_offzet = min (indentation_Ds) + (displacement_peak_to_peak/2);

force_peak_to_peak = max (indentation_Fs)- min (indentation_Fs);
first_force= indentation_Fs(1);
forced_offzet = min (indentation_Fs) + (force_peak_to_peak/2);

displaced_forces = tsmovavg(indentation_Fs,'t',FT_movingaverage,1)-forced_offzet;
u=1;
for ind =1:length(displaced_forces)-1
  if((displaced_forces(ind)>=0) & (displaced_forces(ind+1)<0))
      fzeros_index(u) = ind;
      u=u+1;
  end
end
displaced_displacements = tsmovavg(indentation_Ds,'t',FT_movingaverage,1)-displaced_offzet;
u=1;
for ind =1:length(displaced_displacements)-1
  if((displaced_displacements(ind)>=0) & (displaced_displacements(ind+1)<0))
      dzeros_index(u) = ind;
      u=u+1;
  end
end

for i =1:length(dzeros_index)
    time_shifts(i) = indentation_ts(dzeros_index(i))-indentation_ts(fzeros_index(i));    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%calculating the frequency
frequency_Hz =  (length(fzeros_index)-1)/(indentation_ts(fzeros_index(length(fzeros_index)))-indentation_ts(fzeros_index(1)))
omega= frequency_Hz*2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(g, 'Position'); % grab original on screen size
set(g, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(g,'PaperPositionMode','auto') %set paper pos for printing

%g is for figure 2

%figure(2)

plot(indentation_ts*1000,displaced_forces,'-','LineWidth',1);
set(gca,'fontsize',20)
xlabel('Time (ms)' ,'FontSize', 20);
ylabel('Displaced Value','FontSize', 20);
str101=sprintf('Displaced Values of Force (Blue) and Displacement (Red), Average Time Shift = %f ms, f = %f Hz',mean(time_shifts)*1000,frequency_Hz);
title(str101,'FontSize', 20);
hold on
plot(indentation_ts*1000,displaced_displacements,'-','LineWidth',1);
set(gca,'fontsize',20)
grid on
grid minor

%saveFigure(figure(2),'Displaced_Values');

file_name =strcat(path,'\',Ex_name,'_',my_sweep_index,'_Displaced_Waves')

print(g,file_name,'-dpng', '-r300')
print(g,file_name,'-dtiff')
saveas(g, file_name) % save figure
set(g,'Position', origSize) %set back to original dimensions

to_print_displaced_waves=strcat(file_name,'.png');

%Displaced Values Figure

%figure(3)
% h is for figure(3)
h = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(h, 'Position'); % grab original on screen size
set(h, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(h,'PaperPositionMode','auto') %set paper pos for printing

subplot(2,2,1)
b = plot(Indentation_depth(index:2048),tsmovavg(indentation_Fs,'t',FD_movingaverage,1),'-');
set(gca,'fontsize',20)
xlabel('Indentation Depth (µm)','FontSize', 20);
ylabel('Force (N)','FontSize', 20);
grid on
grid minor    
hold on

subplot(2,2,2);
c=plot(indentation_ts,tsmovavg(indentation_Fs,'t',FT_movingaverage,1),'-');
set(gca,'fontsize',20)
xlabel('Time (s)','FontSize', 20);
ylabel('Force (N)','FontSize', 20);
grid on
grid minor
hold on

subplot(2,2,4);
plot(indentation_ts,Indentation_depth(index:2048),'-');
set(gca,'fontsize',20)

xlabel('Time (s)','FontSize', 20);
ylabel('Indentation Depth (µm)','FontSize', 20); 
hold on

grid minor
grid on

%saveFigure(figure(3),'Indentation_Region');
file_name =strcat(path,'\',Ex_name,'_',my_sweep_index,'_Displaced_Values')

print(h,file_name,'-dpng', '-r300')
print(h,file_name,'-dtiff')
saveas(h, file_name) % save figure
set(h,'Position', origSize) %set back to original dimensions

to_print_displaced_values=strcat(file_name,'.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating The Storage And The Loss Modulus
for i=1:2048
    contact_radius(i) = (Indentation_depth(i)*(2* probe_radius - Indentation_depth(i)))^0.5;
end


for i=1:2048
    contact_area(i) = (probe_radius*Indentation_depth(i)*2*pi);
end

for i =1:2048
    %Storage_Modulus(i) = ((max(FTs)/max(Indentation_depth))* (contact_area(i))^-0.5 *10^12);
    %Loss_Modulus(i) =    ((max(FTs)/max(Indentation_depth))* (contact_area(i))^-0.5 *10^12);
    
    Storage_Modulus(i) =(1-(poisson^2))* ((max(FTs(:,2))/max(Indentation_depth))* (contact_radius(i))^-1 *10^12);
    Loss_Modulus(i)    =(1-(poisson^2))* ((max(FTs(:,2))/max(Indentation_depth))* (contact_radius(i))^-1 *10^12);
end
delta = omega*mean(time_shifts);   
for i =1:2048
    Storagy(i) = Storage_Modulus(i)*cos(delta);
    Lossy(i) =   Loss_Modulus(i)*sin(delta);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k is for figure(4)
k = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(k, 'Position'); % grab original on screen size
set(k, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(k,'PaperPositionMode','auto') %set paper pos for printing


subplot(2,1,1);
plot(Indentation_depth(index:2048),Storagy(index:2048)*10^-3,'LineWidth',2);
set(gca,'fontsize',20)
str1=sprintf('\\omega = %f rad/s, \\delta = %f rad ',omega,delta);
title (str1,'FontSize', 20);
xlabel('Indentation Depth (µm)','FontSize', 20);
ylabel('Storage Modulus (KPa)','FontSize', 20);
grid on
grid minor 
hold on
subplot(2,1,2);
plot(Indentation_depth(index:2048),Lossy(index:2048)*10^-3,'LineWidth',2);
set(gca,'fontsize',20)
xlabel('Indentation Depth (µm)','FontSize', 20);
ylabel('Loss Modulus (KPa)','FontSize', 20);
grid on 
grid minor 

%saveFigure(figure(4),'Storage_Loss_Modulus');
file_name =strcat(path,'\',Ex_name,'_',my_sweep_index,'_Storage_Loss_Modulus')

print(k,file_name,'-dpng', '-r300')
print(k,file_name,'-dtiff')
saveas(k, file_name) % save figure
set(k,'Position', origSize) %set back to original dimensions

to_print_storage_loss=strcat(file_name,'.png');

%%%%%%%saving into file section

cd(path)

mkdir('DMA_Analysis')
new_path =strcat(path,'\','DMA_Analysis');
cd(new_path);

final_analysis_filename= 'DMA_Analysis.csv';
a=exist(final_analysis_filename, 'file')

if(a==0)
	fid = fopen(final_analysis_filename, 'a');
	fprintf(fid, 'Frequency (Hz) ; Amplitude (N); Average Time Shift (ms); Storage Modulus (KPa) ; Loss Modulus (KPa) ; Poisson Ratio');
	fclose(fid)
end


fid = fopen(final_analysis_filename, 'a');
fprintf(fid,strcat('\n',strrep(num2str(frequency_Hz),'.', ','),';',strrep(num2str(max(FDs([1:2048],2))),'.', ','),';',strrep(num2str(mean(time_shifts)*1000),'.', ','),';',strrep(num2str(mean(Storagy(index:2048)*10^-3)),'.', ',') ,';', strrep(num2str(mean(Lossy(index:2048)*10^-3)) ,'.', ','),';', strrep(num2str(poisson) ,'.', ',')));
fclose(fid);



%%%%%%%%%%%%%% opening the file again and plotting the storage/loss modulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dd = dir('*.csv');                           

fileNames = {dd.name}; 
copyfile(dd.name,'temp.csv')
for i=1:length(fileNames)

	file    = memmapfile( 'temp.csv', 'writable', true );
	comma   = uint8(',');
	point   = uint8('.');
	file.Data( transpose( file.Data==comma) ) = point;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plotting and saving storage and loss moduli

data = cell(numel(fileNames),2);
data(:,1) = regexprep(fileNames, '.csv','');

for ii = 1:numel(fileNames)    
   data{ii,2} = dlmread('temp.csv',';',1,0);
   [n,s,r] = xlsread('temp.csv');
end

f = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(f, 'Position'); % grab original on screen size
set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(f,'PaperPositionMode','auto') %set paper pos for printing


frequencies = data{1,2}(:,1)
Storage_Moduli = data{1,2}(:,4);
Loss_Moduli = data{1,2}(:,5);
subplot(2,1,1)

plot(frequencies,Storage_Moduli,'-o','LineWidth',2,'Markersize',10);
grid on;
grid minor;
set(gca,'fontsize',20)
xlabel('Frequency (Hz)','FontSize', 20);
ylabel('Storage Modulus (KPa)','FontSize', 20);
subplot(2,1,2)

plot(frequencies,Loss_Moduli,'-o','LineWidth',2,'Markersize',10);
grid on;
grid minor;
set(gca,'fontsize',20)
xlabel('Frequency (Hz)','FontSize', 20);
ylabel('Loss Modulus (KPa)','FontSize', 20);

file_name =strcat(path,'\',Ex_name,'_Storage_Loss_Moduli')
print(f,file_name,'-dpng', '-r300')
print(f,file_name,'-dtiff')
saveas(f, file_name) % save figure


to_print_storage_loss_moduli=strcat(file_name,'.png');

