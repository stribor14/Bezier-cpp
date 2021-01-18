#ifndef CUSTOMSCENE_H
#define CUSTOMSCENE_H

#include <QGraphicsScene>

#include "qcurve.h"
#include "qpolycurve.h"
#include "qgraphicsviewzoom.h"

namespace Ui
{
class MainWindow;
}

inline QTextStream& qStdOut()
{
  static QTextStream ts(stdout);
  return ts;
}

class CustomScene : public QGraphicsScene
{
private:
  QGraphicsEllipseItem* dot;
  QMap<QGraphicsItem*, QGraphicsLineItem*> line;
  QMap<QGraphicsItem*, QGraphicsLineItem*> tan;
  QMap<QGraphicsItem*, QGraphicsEllipseItem*> byLength;
  bool draw_box_ = false;
  bool draw_inter_ = false;
  bool show_projection = false;
  bool update_curvature = false;
  std::pair<qCurve*, double> t_to_update;
  bool update_cp = false;
  std::pair<QGraphicsItem*, unsigned> cp_to_update;

protected:
  void drawForeground(QPainter* painter, const QRectF& rect) Q_DECL_OVERRIDE;
  void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;
  void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE { Q_UNUSED(mouseEvent) }

  void keyPressEvent(QKeyEvent* keyEvent) Q_DECL_OVERRIDE;
};

#endif // CUSTOMSCENE_H
