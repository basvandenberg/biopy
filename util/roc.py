import numpy
from matplotlib import pyplot
#from sklearn.metrics import roc_curve, auc

'''
Created on Fri Oct 1, 2010

@author: Bastiaan van den Berg
'''

class RocCvCollection(object):

    def __init__(self):
        self.roc_list = []

    def add(self, r):
        self.roc_list.append(r)

    def save_roc_plot(self, f, color='#3465a4', linestyle='-',
                subcolor='#babdb6', sublinestyle='-'):
        fig = self.get_roc_plot(color, linestyle, subcolor, sublinestyle)
        fig.savefig(f)
        pyplot.close(fig)

    def get_roc_plot(self, color='#3465a4', linestyle='-', subcolor='#babdb6',
                sublinestyle='-'):

        # create figure axes
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)

        # plot the random classification line
        ax.plot([0, 1], [0, 1], color='#d3d7cf', linestyle='--')

        # plot all rocs in the collection
        for index, r in enumerate(self.roc_list):
            r.add_to_roc_axes(ax, color=subcolor, linestyle=sublinestyle)

        # plot the average roc
        x, y = self.avg_roc()
        avg, std = self.avg_auc()
        ax.plot(x, y, color=color, label='avg-auc = %0.2f (std = %0.2f))' %
                (avg, std))

        # general plot settings
        ax.grid()
        ax.set_xlabel('false positive rate')
        ax.set_ylabel('true positive rate')
        ax.legend(loc="lower right", prop={'size':8})
        return fig

    def avg_roc(self):
        # TODO implement this, now returning the first ROC as stub
        return self.roc_list[0].roc

    def avg_auc(self):
        '''
        Returns mean and std of the area under the curves of all ROC-curves in
        the collection.

        NOTE: the area under the average roc-curve might be slightly different.
        '''
        aucs = [r.auc() for r in self.roc_list]
        return (numpy.mean(aucs), numpy.std(aucs))

class ROC(object):
    '''
    '''

    def __init__(self, class_labels, predictions, class0=-1, class1=1):
        '''
        ROC-curve for 2-class classification
        The predictions and the class labels must be a list with numbers, 
        negatives for one class, positives for the other class.
        predictions: list with real numbers, negative values for one class
                     positive values for the other
        class_labels: list with class labels, -1.0 for one class, 1.0 for
                      the other class. both values must be present in the
                      list at least once. A ValueError will be raised 
                      otherwise.
        A value error will be raised when the lists have different sizes
        '''
        self._check_args(class_labels, predictions, class0, class1)
        self.predictions = predictions
        self.class_labels = class_labels
        self.class0 = class0
        self.class1 = class1
        self.roc = self._roc()

    def _check_args(self, class_labels, predictions, class0, class1):
        if not(len(predictions) == len(class_labels)):
            raise ValueError('Unequal number of predictions and class labels.')
        # TODO check classes
        #class_set = set(class_labels)
        #if not(val_set == set(classes)):
        #    raise ValueError('Class labels may only be -1.0 and 1.0,' + 
        #                     'and both values must be at least once present' + 
        #                     'in the class labels list.')

    def _roc(self):
        '''
        >>> from bio import roc
        >>> pred = [-1.0, 1.0]
        >>> clab = [-1.0, 1.0]
        >>> r = roc.ROC(pred, clab)
        >>> r.roc
        ([0.0, 0.0, 1.0], [0.0, 1.0, 1.0])

        '''

        # list with x and y values roc-curve
        x = []
        y = []

        class_counts = numpy.bincount(self.class_labels)

        # determine N and P (number of negative and positive class labels)
        n = class_counts[self.class0]
        p = class_counts[self.class1]

        assert(n + p == len(self.class_labels))

        pred_sorted = sorted(self.predictions[:])

        first = pred_sorted[0] - 0.1
        last = pred_sorted[-1] + 0.1
        middle = []
        for i in range(len(pred_sorted) - 1):
            if not(pred_sorted[i] == pred_sorted[i + 1]):
                midpoint = pred_sorted[i] + 0.5 *\
                        (pred_sorted[i + 1] - pred_sorted[i])
                middle.append(midpoint)
        
        thresholds = []
        thresholds.append(first)
        thresholds.extend(middle)
        thresholds.append(last)
        thresholds.reverse()
        
        for threshold in thresholds:

            # determine number of false and true positives
            fp = 0
            tp = 0
            for i in range(len(self.predictions)):
                if(self.predictions[i] > threshold):
                    if(self.class_labels[i] == self.class0):
                        fp += 1
                    elif(self.class_labels[i] == self.class1):
                        tp += 1
                    else:
                        raise ValueError('Labels can only be %i or %i.' %
                                (self.class0, self.class1))

            # calculate false and true positive rate
            fpr = float(fp) / n
            tpr = float(tp) / p
            # add to list
            x.append(fpr)
            y.append(tpr)

        return (x, y)

    def save_roc_plot(self, f, label='', color='#3465a4', linestyle='-'):
        fig = self.get_roc_plot(f, label, color, linestyle)
        fig.savefig(f)
        pyplot.close(fig)

    def get_roc_plot(self, f, label='', color='#3465a4', linestyle='-'):
        fig = pyplot.figure()
        ax = fig.add_subplot(1, 1, 1)
        self.add_to_roc_axes(ax, label, color, linestyle)
        ax.grid()
        ax.set_xlabel('false positive rate')
        ax.set_ylabel('true positive rate')
        ax.legend(loc="lower right")
        return fig

    def add_to_roc_axes(self, ax, label='', color='#babdb6', linestyle='-'):
        x, y = self.roc
        if(len(label) > 0):
            label = '%s (area = %0.2f)' % (label, self.auc())
            ax.plot(x, y, color=color, label=label)
        else:
            ax.plot(x, y, color=color)

    #def auc_roc(self, limit=1.0):
    def auc(self):

        #if(limit < 0.0 or limit > 1.0):
        #    raise ValueError('Limit must be in range [0.0, 1.0].')

        (xvalues, yvalues) = self.roc

        '''
        area = 0.0
        roci = 0

        while roci + 1 < len(xvalues) and xvalues[roci + 1] <= limit:
            x0 = xvalues[roci]
            x1 = xvalues[roci + 1]
            y0 = yvalues[roci]
            y1 = yvalues[roci + 1]
            miny = min(y0, y1)
            binarea = (x1 - x0) * miny + 0.5 * (x1 - x0) * abs(y1 - y0)
            area += binarea
            roci += 1

        # add left part of the last bin (x0 to limit)  
        if(roci + 1 < len(xvalues) and xvalues[roci + 1] > limit):
            x0 = xvalues[roci]
            x1 = xvalues[roci + 1]
            y0 = yvalues[roci]
            y1 = yvalues[roci + 1]
            miny = min(y0, y1)
            xdist = limit - x0
            ydist = (xdist / (x1 - x0)) * abs(y1 - y0)
            binarea = xdist * miny + 0.5 * xdist * ydist
            area = area + binarea

        return area
        '''
        return numpy.trapz(yvalues, xvalues)
