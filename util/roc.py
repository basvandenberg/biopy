import numpy
from matplotlib import pyplot
from sklearn.metrics import roc_curve, auc

'''
Created on Fri Oct 1, 2010

@author: Bastiaan van den Berg
'''

class ROC(object):
    '''
    '''

    def __init__(self, predictions, class_labels):
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
        self._check_args(predictions, class_labels)
        self.predictions = predictions
        self.class_labels = class_labels  
        self.roc = self._roc()

    def _check_args(self, predictions, class_labels):
        if not(len(predictions) == len(class_labels)):
            raise ValueError('Unequal number of predictions and class labels.')
        val_set = set(class_labels)
        if not(val_set == set([-1.0, 1.0])):
            raise ValueError('Class labels may only be -1.0 and 1.0,' + 
                             'and both values must be at least once present' + 
                             'in the class labels list.')

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

        # determine N and P (number of negative and positive class labels)
        n = self.class_labels.count(-1.0)
        p = self.class_labels.count(1.0)

        assert(n + p == len(self.class_labels))

        pred_sorted = sorted(self.predictions[:])

        first = pred_sorted[0] - 0.1
        last = pred_sorted[-1] + 0.1
        middle = []
        for i in range(len(pred_sorted) - 1):
            if not(pred_sorted[i] == pred_sorted[i + 1]):
                midpoint = pred_sorted[i] + 0.5 * (pred_sorted[i + 1] - pred_sorted[i])
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
                    if(self.class_labels[i] == -1.0):
                        fp += 1
                    elif(self.class_labels[i] == 1.0):
                        tp += 1
                    else:
                        raise ValueError('Labels can only be -1.0 or 1.0.')

            # calculate false and true positive rate
            fpr = float(fp) / n
            tpr = float(tp) / p
            # add to list
            x.append(fpr)
            y.append(tpr)

        return (x, y)

    def plot_roc(self, f):
        (x, y) = self.roc
        fig = pyplot.figure()
        pyplot.plot(x, y)
        pyplot.grid(True)
        fig.savefig(f)
        fig.clear()

    def auc_roc(self, limit=1.0):

        if(limit < 0.0 or limit > 1.0):
            raise ValueError('Limit must be in range [0.0, 1.0].')

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
