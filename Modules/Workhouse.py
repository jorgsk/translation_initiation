import rpy
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

rpy.r.options(warn=-1) # suppress warning messages from R (spearman equal ranks)

def ReadingFrame(x):
    """ Assume the sequences begin with reading frames 0, 1, and 2, and give
    which reading frame x belongs to. """
    frames = [(x-n)%3 for n in range(3)]
    frame = frames.index(0)
    return frame

def SubSeqLocater(string, keyword):
    """ Returns the (start,stop) indices of all the subsequences of keyword in
    string. Start and stop are given as in slice notation (from and including
    start to but not including stop -- [start, stop)). """

    locations = []
    strlen = len(string)
    keylen = len(keyword)
    last   = 0
    while last < strlen:
        start = string.find(keyword, last)
        if start == -1: break
        stop  = start + keylen
        last  = stop
        locations.append((start,stop))
    return locations

def StringOrFloat(incoming):
    """ Return float if element is float, otherwise return unmodified. """
    datatype = type(incoming)
    if datatype == str:
        try:
            fl00t = float(incoming)
            return fl00t
        except:
            return incoming
    elif (datatype is list) or (datatype is np.ndarray):
        outgoing = []
        for element in incoming:
            try:
                element = float(element)
                outgoin.append(element)
            except:
                outgoing.append(element)
        return outgoing
    raise ValueError('Input must be string/float/int or list/array of these. ')

def NumberFormatter(indata, sign):
    """ Accepts floats/ints or arrays of floats/ints. Converts these to strings
    of with as many digits as given in sign. If in list/array, run program
    RECURSIVELY(!) to return list/array. Strings are returned as they come in. """
    # '{0:.2g}'.format(0.0007478792348) #Out[40]: '0.00075'
    datatype = type(indata)

    if datatype in [int, np.int32, np.int64, float, np.float64, np.float32]:
        formstr = '{0:.'+str(sign)+'g}'
        return formstr.format(indata)

    if datatype in [list, np.ndarray]:
        for index, value in enumerate(indata[:]):
            indata[index] = NumberFormatter(value,sign)
        return indata

    if datatype is str:
        return indata

    raise ValueError('Input must be int/float/str or list/array of these. ')

def DataSetCorrelations(dataset):
    """ Return correlations between all lists in 'dataset'."""
    result = []
    for indx1, row1 in enumerate(dataset):
        for indx2, row2 in enumerate(dataset):
            (corr, pval) = Spearman(row1,row2)
            result.append([corr, pval, (indx1,indx2)])
    return result

def Spearman(x, y):
    """ Return Spearman (rho, pval) from lists x and y """
    test = rpy.r.cor_test(x, y, method="spearman")
    rho = test['estimate']['rho']
    pval = test['p.value']
    return rho, pval

def DictSort(d, xx, rev=True):
    """ Return dictionary d sorted by keys (xx=0) or values (xx=1). Reverse
    order (largest first) is default """
    return sorted(d.iteritems(), key=itemgetter(xx), reverse=rev)

def StdPlotter(xdata,ydata,ystdv,xlab='xlabel',ylab='ylabel', title='title'):
    """ Scatter plot with error bars on the y-data. """
    fig = plt.figure()
    ax = fig.add_subplot(111) # to have a plot with axes to work on! not just a figure.
    ax.errorbar(xdata,ydata,yerr=ystdv,fmt='ro')
    ax.set_xlabel(xlab, fontsize=20)
    ax.set_ylabel(ylab, fontsize=20)
    ax.set_title(title, fontsize=20)
    ax.set_ylim(-2,12)
    plt.show()

def RelativeStandardDeviation(data,data_stdv):
    """ Calculate the relative standard deviations stdv of the set py """
    data = np.array(data)
    stdv = np.array(data_stdv)
    relative = np.zeros(43)
    for i in range(len(data)):
        temp = 1-(data[i]-stdv[i])/data[i]
        relative[i] = temp
    plt.scatter(data,relative)
    plt.xlabel('Productive yield', fontsize=15)
    plt.ylabel('Relative standard deviation', fontsize=15)
    plt.title(r'Relative std: $\sigma_{\text{rel}} = \frac{\mu -'
              '\sigma}{\mu}$. Rank correlation: R = -0.17, p = 0.26', fontsize=15)
    spearStats = Spearman(data,relative)
    print('Correlation coefficient: '+str(spearStats[0]),'p-value: '+str(speaStats[1]))

