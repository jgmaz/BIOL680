fname = 'R016-2012-10-08-CSC03b.Ncs';
csc_new = LoadCSCnew(fname);
csc_info_new = getHeader(csc_new);
csc_old = LoadCSC(fname);

subplot(221)
plot(csc_new)
title('New LoadCSC')
subplot(222)
plot(csc_old)
title('Old LoadCSC')
subplot(223)
plot(diff(Range(csc_new)),'r')
subplot(224)
plot(diff(Range(csc_old)),'r')